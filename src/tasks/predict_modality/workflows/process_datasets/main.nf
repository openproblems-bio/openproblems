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
      key: "check_dataset_schema_rna",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input_rna")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input_rna,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset_rna": checks["exit_code"] == 0 ? state.input_rna : null,
        ]
      }
    )

    | check_dataset_schema.run(
      key: "check_dataset_schema_other_mod",
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input_other_mod")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input_other_mod,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset_other_mod": checks["exit_code"] == 0 ? state.input_other_mod : null,
        ]
      }
    )
    | view{"test: ${it}"}

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
