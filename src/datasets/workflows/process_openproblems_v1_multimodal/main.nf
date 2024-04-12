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
      args: [normalization_id: "log_cp10k", n_cp: 10000]
    ),
    log_cp.run(
      key: "log_cpm",
      args: [normalization_id: "log_cpm", n_cp: 1000000]
    ),
    sqrt_cp.run(
      key: "sqrt_cp10k",
      args: [normalization_id: "sqrt_cp10k", n_cp: 10000]
    ),
    sqrt_cp.run(
      key: "sqrt_cpm",
      args: [normalization_id: "sqrt_cpm", n_cp: 1000000]
    ),
    l1_sqrt.run(
      key: "l1_sqrt",
      args: [normalization_id: "l1_sqrt"]
    ),
    log_scran_pooling.run(
      key: "log_scran_pooling",
      args: [normalization_id: "log_scran_pooling"]
    )
  ]

  output_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // fetch data from legacy openproblems
    | openproblems_v1_multimodal.run(
      fromState: [
        "input_id": "input_id",
        "obs_cell_type": "obs_cell_type",
        "obs_batch": "obs_batch",
        "obs_tissue": "obs_tissue",
        "layer_counts": "layer_counts",
        "sparse": "sparse",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description",
        "dataset_organism": "dataset_organism"
      ],
      toState: [
        "raw_mod1": "output_mod1",
        "raw_mod2": "output_mod2"
      ]
    )

    // subsample if need be
    | subsample.run(
      runIf: { id, state -> state.do_subsample },
      fromState: [
        "input": "raw_mod1",
        "input_mod2": "raw_mod2",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_cell_type_categories": "keep_cell_type_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      toState: [
        "raw_mod1": "output",
        "raw_mod2": "output_mod2"
      ]
    )

    // run normalization methods
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
      fromState: ["input": "raw_mod1"],
      toState: { id, output, state, comp -> 
        state + [
          "normalization_id": comp.name,
          "normalized_mod1": output.output
        ]
      }
    )

    // run normalization methods on second modality
    // TODO: can we change this to DSB?
    | prot_clr.run(
      runIf: { id, state -> state.mod2 == "ADT" },
      args: [normalization_id: "prot_clr"],
      fromState: ["input": "raw_mod2"],
      toState: ["normalized_mod2": "output"]
    )
    | atac_tfidf.run(
      runIf: { id, state -> state.mod2 == "ATAC" },
      args: [normalization_id: "atac_tfidf"],
      fromState: ["input": "raw_mod2"],
      toState: ["normalized_mod2": "output"]
    )

    | svd.run(
      fromState: [
        "input": "normalized_mod1",
        "input_mod2": "normalized_mod2"
      ],
      toState: [
        "svd_mod1": "output",
        "svd_mod2": "output_mod2"
      ]
    )

    | hvg.run(
      fromState: [ "input": "svd_mod1" ],
      toState: [ "hvg_mod1": "output" ]
    )

    | hvg.run(
      fromState: [ "input": "svd_mod2" ],
      toState: [ "hvg_mod2": "output" ]
    )

    // add synonyms
    | map{ id, state ->
      [id, state + [
        "output_dataset_mod1": state.hvg_mod1,
        "output_dataset_mod2": state.hvg_mod2
      ]]
    }

    | extract_metadata.run(
      key: "extract_metadata_mod1",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "output_dataset_mod1")
        // workaround: convert GString to String
        schema = iterateMap(schema, { it instanceof GString ? it.toString() : it })
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.output_dataset_mod1,
          "schema": schemaYaml
        ]
      },
      toState: ["output_meta_mod1": "output"]
    )

    | extract_metadata.run(
      key: "extract_metadata_mod2",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "output_dataset_mod2")
        // workaround: convert GString to String
        schema = iterateMap(schema, { it instanceof GString ? it.toString() : it })
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.output_dataset_mod2,
          "schema": schemaYaml
        ]
      },
      toState: ["output_meta_mod2": "output"]
    )
    
    // only output the files for which an output file was specified
    | setState([
      "output_dataset_mod1",
      "output_dataset_mod2",
      "output_meta_mod1",
      "output_meta_mod2",
      "_meta"
    ])

  emit:
  output_ch
}
