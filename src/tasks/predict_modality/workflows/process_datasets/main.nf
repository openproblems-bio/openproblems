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
  output_ch = input_ch

    | check_dataset_schema.run(
      key: "check_dataset_schema_mod1",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input_mod1")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input_mod1,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset_mod1": checks["exit_code"] == 0 ? state.input_mod1 : null,
        ]
      }
    )

    | check_dataset_schema.run(
      key: "check_dataset_schema_mod2",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input_mod2")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input_mod2,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset_mod2": checks["exit_code"] == 0 ? state.input_mod2 : null,
        ]
      }
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset_mod1 != null &&
      state.dataset_mod2 != null
    }

    | process_dataset.run(
      fromState: [
        input_mod1: "dataset_mod1",
        input_mod2: "dataset_mod2",
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
