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
  output_ch = input_ch

    // TODO: check schema based on the values in `config`
    // instead of having to provide a separate schema file
    | check_dataset_schema.run(
      key: "check_dataset_schema_rna",
            fromState: { id, state ->
        // as a resource
        [
          "input": state.input_rna,
          "schema": meta.resources_dir.resolve("file_common_dataset_rna.yaml")
        ]
      },
      args: [
        "stop_on_error": false
      ],
      toState: [
        "dataset_rna": "output",
        "dataset_checks": "checks"
      ]
    )

    | check_dataset_schema.run(
      key: "check_dataset_schema_other_mod",
            fromState: { id, state ->
        // as a resource
        [
          "input": state.input_other_mod,
          "schema": meta.resources_dir.resolve("file_common_dataset_other_mod.yaml")
        ]
      },
      args: [
        "stop_on_error": false
      ],
      toState: [
        "dataset_other_mod": "output",
        "dataset_checks": "checks"
      ]
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset_rna != null &&
      state.dataset_other_mod != null
    }

    | process_dataset.run(
      fromState: [
        input_rna: "dataset_rna",
        input_other_mod: "dataset_other_mod",
        output_train_mod1: "output_train_mod1",
        output_train_mod2: "output_train_mod2",
        output_test_mod1: "output_test_mod1",
        output_test_mod2: "output_test_mod2"
      ],
      toState: [
        "output_train_mod1",
        "output_train_mod2",
        "output_test_mod1",
        "output_test_mod2"
      ]
    )

    // only output the files for which an output file was specified
    | setState ([
        "output_train_mod1",
        "output_train_mod2",
        "output_test_mod1",
        "output_test_mod2"
      ])

  emit:
  output_ch
}
