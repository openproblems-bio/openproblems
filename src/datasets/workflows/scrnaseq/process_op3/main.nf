include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  // create different normalization methods by overriding the defaults
  normalization_methods = [
    log_cp.run(
      key: "log_cp10k",
      args: [normalization_id: "log_cp10k", n_cp: 10000],
    ),
    log_cp.run(
      key: "log_cpm",
      args: [normalization_id: "log_cpm", n_cp: 1000000],
    ),
    sqrt_cp.run(
      key: "sqrt_cp10k",
      args: [normalization_id: "sqrt_cp10k", n_cp: 10000],
    ),
    sqrt_cp.run(
      key: "sqrt_cpm",
      args: [normalization_id: "sqrt_cpm", n_cp: 1000000],
    ),
    l1_sqrt.run(
      key: "l1_sqrt",
      args: [normalization_id: "l1_sqrt"],
    ),
    log_scran_pooling.run(
      key: "log_scran_pooling",
      args: [normalization_id: "log_scran_pooling"],
    )
  ]

  output_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // fetch data from OP3 dataset
    | openproblems_op3.run(
      fromState: [
	"input": "input",
        "data_type": "data_type",
        "donor_id": "donor_id",
        "cell_type": "cell_type",
        "perturbation": "perturbation",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description"
      ],
      toState: ["output_raw": "output"]
    )

    // subsample if so desired
    | subsample.run(
      runIf: { id, state -> state.do_subsample },
      fromState: [
        "input": "output_raw",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_cell_type_categories": "keep_cell_type_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      args: [output_mod2: null],
      toState: ["output_raw": "output"]
    )

    | runEach(
      components: normalization_methods,
      id: { id, state, comp -> 
        if (state.normalization_methods.size() > 1) {
          id + "/" + comp.name
        } else {
          id
        }
      },
      filter: { id, state, comp ->
        comp.name in state.normalization_methods
      },
      fromState: ["input": "output_raw"],
      toState: { id, output, state, comp ->
        state + [
          output_normalized: output.output,
          normalization_id: comp.name
        ]
      }
    )

    | hvg.run(
      fromState: ["input": "output_normalized"],
      toState: ["output_hvg": "output"]
    )

    | pca.run(
      fromState: ["input": "output_hvg"],
      toState: ["output_pca": "output" ]
    )

    | knn.run(
      fromState: ["input": "output_pca", "output_compression": "output_compression"],
      toState: ["output_knn": "output"]
    )

    // add synonym
    | map{ id, state ->
      [id, state + [output_dataset: state.output_knn]]
    }

    | extract_uns_metadata.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "output_dataset")
        // workaround: convert GString to String
        schema = iterateMap(schema, { it instanceof GString ? it.toString() : it })
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.output_dataset,
          "schema": schemaYaml
        ]
      },
      toState: ["output_meta": "output"]
    )

    // only output the files for which an output file was specified
    | setState([
      "output_dataset",
      "output_meta",
      "output_raw",
      "output_normalized",
      "output_pca",
      "output_hvg",
      "output_knn",
      "_meta"
    ])

  emit:
  output_ch
} 
