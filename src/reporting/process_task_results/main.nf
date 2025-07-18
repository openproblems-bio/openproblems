workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | get_task_info.run(
      fromState: [
        "input": "input_task_info"
      ],
      toState: ["output_task": "output"]
    )

    // extract task id from task info
    | map { id, state ->
      def task_id = readJson(state.output_task).task_id
      [id, state + ["task_id": task_id]]
    }

    | get_dataset_info.run(
      fromState: [
        "input": "input_dataset_info",
      ],
      toState: ["output_dataset": "output"]
    )

    | get_method_info.run(
      fromState: [
        "input": "input_method_configs",
      ],
      toState: ["output_method": "output"]
    )

    | get_metric_info.run(
      fromState: [
        "input": "input_metric_configs",
      ],
      toState: ["output_metric": "output"]
    )

    | get_results.run(
      fromState: [
        "input_scores": "input_scores",
        "input_trace": "input_trace",
        "input_dataset_info": "output_dataset",
        "input_method_info": "output_method",
        "input_metric_info": "output_metric"
      ],
      toState: [
        "output_results": "output"
      ]
    )

    | generate_qc.run(
      fromState: [
        "input_task_info": "output_task",
        "input_dataset_info": "output_dataset",
        "input_method_info": "output_method",
        "input_metric_info": "output_metric",
        "input_results": "output_results"
      ],
      toState: ["output_qc": "output"]
    )

    | combine_output.run(
      fromState: [
        "input_task_info": "output_task",
        "input_quality_control": "output_qc",
        "input_metric_info": "output_metric",
        "input_method_info": "output_method",
        "input_dataset_info": "output_dataset",
        "input_results": "output_results"
      ],
      toState: ["output_combined": "output"]
    )

    | setState([
      "output_combined": "output_combined",
      "output_task_info": "output_task",
      "output_dataset_info": "output_dataset",
      "output_method_info": "output_method",
      "output_metric_info": "output_metric",
      "output_results": "output_results",
      "output_quality_control": "output_qc"
    ])

  emit:
  output_ch
}
