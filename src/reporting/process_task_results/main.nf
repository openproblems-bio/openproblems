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
        "input_metric_info": "output_metric",
        "datasets_include": "datasets_include",
        "datasets_exclude": "datasets_exclude",
        "methods_include": "methods_include",
        "methods_exclude": "methods_exclude",
        "metrics_include": "metrics_include",
        "metrics_exclude": "metrics_exclude"
      ],
      toState: [
        "output_results": "output"
      ]
    )

    | filter_results.run(
      runIf: { id, state ->
        // Only run filtering if there are include/exclude lists defined
        return state.datasets_exclude || state.methods_exclude || state.metrics_exclude ||
          state.datasets_include || state.methods_include || state.metrics_include
      },
      fromState: [
        "input_dataset_info": "output_dataset",
        "input_method_info": "output_method",
        "input_metric_info": "output_metric",
        "input_results": "output_results",
        "datasets_include": "datasets_include",
        "datasets_exclude": "datasets_exclude",
        "methods_include": "methods_include",
        "methods_exclude": "methods_exclude",
        "metrics_include": "metrics_include",
        "metrics_exclude": "metrics_exclude"
      ],
      toState: [
        "output_dataset": "output_dataset_info",
        "output_method": "output_method_info",
        "output_metric": "output_metric_info",
        "output_results": "output_results"
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

    | render_report.run(
      fromState: [
        "input_task_results": "output_combined"
      ],
      toState: ["output_report": "output"]
    )

    | setState([
      "output_data": "output_combined",
      "output_report": "output_report"
    ])

  emit:
  output_ch
}
