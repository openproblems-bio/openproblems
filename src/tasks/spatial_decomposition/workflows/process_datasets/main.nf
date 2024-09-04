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
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset": checks["exit_code"] == 0 ? state.input : null,
        ]
      }
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset != null
    }

    | dataset_simulator.run(
      runIf: {id, state -> state.alpha}, 
      fromState: [ 
        input: "dataset", 
        alpha: "alpha"
      ],
      toState: [ dataset: "simulated_data"], 
      auto: [publish: true]
    )

    | process_dataset.run(
      fromState: [ input: "dataset" ],
      toState: [
        output_single_cell: "output_single_cell",
        output_spatial_masked: "output_spatial_masked",
        output_solution: "output_solution" 
      ]
    )

    // only output the files for which an output file was specified
    | setState(["output_single_cell", "output_spatial_masked", "output_solution"])

  emit:
  output_ch
}
