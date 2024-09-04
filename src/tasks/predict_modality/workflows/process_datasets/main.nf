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
  
    // Check if the input datasets match the desired format --------------------------------
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

    // Use datasets in both directions (mod1 -> mod2 and mod2 -> mod1) ---------------------
    // extract the dataset metadata
    | extract_metadata.run(
      key: "extract_metadata",
      fromState: [input: "dataset_mod1"],
      toState: { id, output, state ->
        def uns = readYaml(output.output).uns
        state + [
          "dataset_id": uns.dataset_id,
          "normalization_id": uns.normalization_id
        ]
      }
    )

    // Add swap direction to the state and set new id
    | flatMap{id, state -> 
      ["normal", "swap"].collect { dir ->
        // Add direction (normal / swap) to id  
        // Note: this id is added before the normalisation id  
        // Example old id: dataset_loader/dataset_id/normalization_id  
        // Example new id: dataset_loader/dataset_id/direction/normalization_id
        def orig_dataset_id = id.replaceAll("/${state.normalization_id}", "")
        def normalization_id = id.replaceAll("^${orig_dataset_id}", "")
        def new_dataset_id = orig_dataset_id + "/" + dir
        def new_id = new_dataset_id + normalization_id

        [new_id, state + [dataset_id: new_dataset_id, direction: dir, "_meta": [join_id: id]]]
      }
    }

    | process_dataset.run(
      fromState: { id, state ->
        def swap_state = state.direction == "swap" ? true : false
        [
          dataset_id: state.dataset_id,
          input_mod1: state.dataset_mod1,
          input_mod2: state.dataset_mod2,
          output_train_mod1: state.output_train_mod1,
          output_train_mod2: state.output_train_mod2,
          output_test_mod1: state.output_test_mod1,
          output_test_mod2: state.output_test_mod2,
          swap: swap_state
        ]
      },
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
      "output_test_mod2",
      "_meta"
    ])

  emit:
  output_ch
}
