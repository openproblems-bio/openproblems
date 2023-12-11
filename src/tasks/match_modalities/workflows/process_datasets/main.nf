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
      key: "check_dataset_schema_mod1",
      fromState: { id, state ->
        // as a resource
        [
          "input": state.input,
          "schema": meta.resources_dir.resolve("file_common_dataset_mod1.yaml")
        ]
      },
      args: [
        "stop_on_error": false
      ],
      toState: [
        "dataset_mod1": "output"
      ]
    )
    | check_dataset_schema.run(
      key: "check_dataset_schema_mod2",
      fromState: { id, state ->
        // as a resource
        [
          "input": state.input,
          "schema": meta.resources_dir.resolve("file_common_dataset_mod2.yaml")
        ]
      },
      args: [
        "stop_on_error": false
      ],
      toState: [
        "dataset_mod2": "output"
      ]
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset_mod1 != null && state.dataset_mod2 != null
    }

    | process_dataset.run(
      fromState: [ input_mod1: "dataset_mod1", input_mod2: "dataset_mod2" ],
      toState: [
        "output_mod1",
        "output_mod2",
        "output_solution_mod1",
        "output_solution_mod2"
      ]
    )

    // only output the files for which an output file was specified
    | setState([
      "output_mod1",
      "output_mod2",
      "output_solution_mod1",
      "output_solution_mod2"
    ])

  emit:
  output_ch
}
