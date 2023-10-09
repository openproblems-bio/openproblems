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
    l1_sqrt,
    log_scran_pooling
  ]

  output_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // fetch data from legacy openproblems
    | openproblems_v1.run(
      fromState: [
        "dataset_id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts",
        "sparse", "dataset_name", "data_url", "data_reference", "dataset_summary",
        "dataset_description", "dataset_organism"
      ],
      toState: ["raw": "output"]
    )
    
    // subsample if so desired
    | subsample.run(
      runIf: { id, state -> state.do_subsample },
      fromState: [
        "input": "raw",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_celltype_categories": "keep_celltype_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      args: [output_mod2: null],
      toState: [raw: "output"]
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
      fromState: ["input": "raw"],
      toState: ["normalized": "output"]
    )

    | pca.run(
      fromState: ["input": "normalized"],
      toState: ["pca": "output" ]
    )

    | hvg.run(
      fromState: ["input": "pca"],
      toState: ["hvg": "output"]
    )

    | knn.run(
      fromState: ["input": "hvg"],
      toState: ["knn": "output"]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.knn,
          checks: null
        ]
      },
      toState: ["dataset": "output", "meta": "meta"]
    )

    // only output the files for which an output file was specified
    | setState{ id, state ->
      [
        "output_dataset": state.output_dataset ? state.dataset : null,
        "output_meta": state.output_meta ? state.meta : null,
        "output_raw": state.output_raw ? state.raw : null,
        "output_normalized": state.output_normalized ? state.normalized : null,
        "output_pca": state.output_pca ? state.pca : null,
        "output_hvg": state.output_hvg ? state.hvg : null,
        "output_knn": state.output_knn ? state.knn : null,
        "_meta": state._meta
      ]
    }

  emit:
  output_ch
}