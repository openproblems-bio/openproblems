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

    // fetch data from legacy openproblems
    | openproblems_v1.run(
      fromState: [
        "dataset_id": "id",
        "obs_celltype": "obs_celltype",
        "obs_batch": "obs_batch",
        "obs_tissue": "obs_tissue",
        "layer_counts": "layer_counts",
        "sparse": "sparse",
        "dataset_name": "dataset_name",
        "data_url": "data_url",
        "data_reference": "data_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description",
        "dataset_organism": "dataset_organism",
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
        "keep_celltype_categories": "keep_celltype_categories",
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
      toState: { id, state, output, comp ->
        state + [
          output_normalization: output.output,
          normalization_id: comp.name
        ]
      }
    )

    | pca.run(
      fromState: ["input": "output_normalized"],
      toState: ["output_pca": "output" ]
    )

    | hvg.run(
      fromState: ["input": "output_pca"],
      toState: ["output_hvg": "output"]
    )

    | knn.run(
      fromState: ["input": "output_hvg"],
      toState: ["output_knn": "output"]
    )

    | check_dataset_schema.run(
      fromState: ["input": "output_knn"],
      toState: ["output_dataset": "output", "output_meta": "meta"]
    )

    | filter{ id, state ->
      def uns = (new org.yaml.snakeyaml.Yaml().load(state.output_meta)).uns
      def expected_id = "${uns.dataset_id}/${uns.normalization_id}"
      
      def is_ok = id == expected_id
      
      if (!is_ok) {
        println("DETECTED ID MISMATCH: $id != $expected_id.\nTuple:\n${toYamlBlob([id, state])}\n")
      }

      is_ok
    }

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