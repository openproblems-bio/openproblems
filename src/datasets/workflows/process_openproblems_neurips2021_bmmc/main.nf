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

    | decompress_gzip.run(
      fromState: ["input": "input"],
      toState: ["input_decompressed": "output"]
    )

    // process neurips downloaded dataset
    | openproblems_neurips2021_bmmc.run(
      fromState: [
        "input": "input_decompressed",
        "mod1": "mod1",
        "mod2": "mod2",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description",
        "dataset_organism": "dataset_organism"
      ],
      toState: [
        "raw_rna": "output_mod1",
        "raw_other_mod": "output_mod2"
      ]
    )

    // subsample if need be
    | subsample.run(
      runIf: { id, state -> state.do_subsample },
      fromState: [
        "input": "raw_rna",
        "input_mod2": "raw_other_mod",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_cell_type_categories": "keep_cell_type_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      toState: [
        "raw_rna": "output",
        "raw_other_mod": "output_mod2"
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
      fromState: ["input": "raw_rna"],
      toState: { id, output, state, comp -> 
        state + [
          "normalization_id": comp.name,
          "normalized_rna": output.output
        ]
      }
    )

    // run normalization methods on second modality
    | log_cp.run(
      key: "log_cp10k_adt",
      runIf: { id, state -> state.mod2 == "ADT" },
      args: [normalization_id: "log_cp10k", n_cp: 10000]
      fromState: ["input": "raw_other_mod"],
      toState: ["normalized_other_mod": "output"]
    )
    | normalization_methods[0].run( // TODO: change this normalization method
      key: "log_cp10k_atac",
      runIf: { id, state -> state.mod2 == "ATAC" },
      args: [normalization_id: "log_cp10k", n_cp: 10000]
      fromState: ["input": "raw_other_mod"],
      toState: ["normalized_other_mod": "output"]
    )

    | svd.run(
      fromState: [
        "input": "normalized_rna",
        "input_mod2": "normalized_other_mod"
      ],
      toState: [
        "svd_rna": "output",
        "svd_other_mod": "output_mod2"
      ]
    )

    | hvg.run(
      fromState: [ "input": "svd_rna" ],
      toState: [ "hvg_rna": "output" ]
    )

    // TODO: should this only run on ATAC? or even not at all?s
    | hvg.run(
      key: "hvg_other_mod",
      fromState: [ "input": "svd_other_mod" ],
      toState: [ "hvg_other_mod": "output" ]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          "input": state.hvg_rna,
          "checks": null
        ]
      },
      toState: [
        "output_rna": "output",
        "output_meta_rna": "meta"
      ]
    )

    | check_dataset_schema.run(
      key: "check_dataset_schema_other_mod",
      fromState: { id, state ->
        [
          "input": state.hvg_other_mod,
          "checks": null
        ]
      },
      toState: [
        "output_other_mod": "output",
        "output_meta_other_mod": "meta"
      ]
    )

    // only output the files for which an output file was specified
    | setState([
      "output_rna",
      "output_other_mod",
      "output_meta_rna",
      "output_meta_other_mod",
      "_meta"
    ])

  emit:
  output_ch
}
