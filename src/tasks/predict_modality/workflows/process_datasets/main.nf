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

  direction = Channel.of("normal", "swap")

  output_ch = input_ch
  
    | combine(direction)

    // Add swap direction to the state and set new id
    | map{id, state, dir -> 
      // Add direction (normal / swap) to id  
      // Note: this id is added before the normalisation id  
      // Example old id: dataset_loader/dataset_id/normalization_id  
      // Example new id: dataset_loader/dataset_id_direction/normalization_id
      def id_split = id.tokenize("/")
      def norm = id_split.takeRight(1)[0]
      def new_id = id_split.dropRight(1).join("/") + "_" + dir + "/" + norm
      
      [new_id, state + [direction: dir, "_meta": [join_id: id]]]
    }

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
      fromState: { id, state ->
        def swap_state = state.direction == "swap" ? true : false
        [
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

    // extract the dataset metadata
    | extract_metadata.run(
      key: "extract_metadata_mod1",
      fromState: [input: "output_train_mod1"],
      toState: { id, output, state ->
        state + [
          modality_mod1: readYaml(output.output).uns.modality
        ]
      }
    )

    // extract the dataset metadata
    | extract_metadata.run(
      key: "extract_metadata_mod2",
      fromState: [input: "output_train_mod2"],
      toState: { id, output, state ->
        state + [
          modality_mod2: readYaml(output.output).uns.modality
        ]
      }
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
