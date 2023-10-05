workflow auto {
  findStates(params, config)
    | run_wf
    | publishStates([:])
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // TODO: check schema based on the values in `config`
    // instead of having to provide a separate schema file
    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.input,
          schema: state.dataset_schema,
          output: '$id.$key.output.h5ad',
          stop_on_error: false,
          checks: null
        ]
      },
      toState: { id, output, state ->
        state + [ dataset: output.output ]
      }
    )

    | filter { id, state ->
      state.dataset != null
    }

    | process_dataset.run(
      fromState: [
        input: "dataset",
        output_dataset: "output_dataset",
        output_solution: "output_solution"
      ],
      toState: [dataset: "output_dataset", solution: "output_solution"]
    )

    // only output the files for which an output file was specified
    | setState { id, state ->
      [
        "output_dataset": state.output_dataset ? state.dataset : null,
        "output_solution": state.output_solution ? state.solution : null
      ]
    }

  emit:
  output_ch
}
