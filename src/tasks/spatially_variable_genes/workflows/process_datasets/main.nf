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

    | select_reference.run(
      fromState: [
        input: "dataset",
        num_features: "num_reference_genes",
        coord_type_proc: "coord_type_proc"
      ],
      toState: [dataset: "output"]
    )

    | simulate_svg.run(
      fromState: [
        input: "dataset",
        gp_k: "gp_k_sim",
        select_top_variable_genes: "select_top_variable_genes_sim"
      ],
      toState: [
        dataset_simulated: "output"
      ]
    )

    | log_cp.run(
      fromState: [
        input: "dataset_simulated",
      ],
      toState: [
        dataset_simulated_normalized: "output"
      ],
      args: [n_cp: -1]
    )

    | split_dataset.run(
      fromState: [
        input: "dataset_simulated_normalized"
      ],
      toState: [
        output_dataset: "output_dataset",
        output_solution: "output_solution" 
      ]
    )

    // only output the files for which an output file was specified
    | setState(["output_dataset", "output_solution", "dataset_simulated_normalized"])

  emit:
  output_ch
}
