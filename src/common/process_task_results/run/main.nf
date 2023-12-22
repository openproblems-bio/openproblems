// workflow auto {
//   findStates(params, meta.config)
//     | meta.workflow.run(
//       auto: [publish: "state"]
//     )
// }

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | get_task_info.run(
      key: "task_info",
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

    | get_method_info.run(
      fromState: [ 
        "input": "input_method_configs",
        "task_id" : "task_id"
      ],
      toState: ["output_method": "output"]
    )

    | get_metric_info.run(
      fromState: [ 
        "input": "input_metric_configs",
        "task_id" : "task_id"
      ],
      toState: ["output_metric": "output"]
    )

    | get_dataset_info.run(
      fromState: [
        "task_id" : "task_id",
        "input": "input_dataset_info",
      ],
      toState: ["output_dataset": "output"]
    )

    | get_results.run(
      fromState: [ 
        "task_id": "task_id",
        "input_scores": "input_scores",
        "input_execution" : "input_execution"
      ],
      toState: ["output_results": "output"]
    )

    | generate_qc.run(
      fromState: [
        "task_info": "output_task",
        "method_info": "output_method",
        "metric_info": "output_metric",
        "dataset_info": "output_dataset",
        "results": "output_results"
      ],
      toState: ["output_qc": "output"]
    )

    | map{ id, state ->
      def new_state = [
        output_scores: state.output_results,
        output_method_info: state.output_method,
        output_metric_info: state.output_metric,
        output_dataset_info: state.output_dataset,
        output_task_info: state.output_task,
        output_qc: state.output_qc
      ]

      [id, new_state]
    }

  emit:
  output_ch
}