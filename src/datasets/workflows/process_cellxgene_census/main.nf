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

    // fetch data from legacy openproblems
    | cellxgene_census.run(
      fromState: [
        "input_uri": "input_uri",
        "census_version": "census_version",
        "species": "species",
        "obs_value_filter": "obs_value_filter",
        "cell_filter_grouping": "cell_filter_grouping",
        "cell_filter_minimum_count": "cell_filter_minimum_count",
        "obs_batch": "obs_batch",
        "obs_batch_separator": "obs_batch_separator",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
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