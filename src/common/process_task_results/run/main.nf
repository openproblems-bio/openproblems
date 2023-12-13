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

    | get_method_info.run(
      fromState: [ 
        "input": "input_method_configs",
        "task_id" : "task_id",
        "output": "output_method_info"
      ],
      toState: { id, output, state ->
        state + [output_method: output.output]
      }
    )

    | get_metric_info.run(
      fromState: [ 
        "input": "input_metric_configs",
        "task_id" : "task_id",
        "output": "output_metric_info"
      ],
      toState: { id, output, state ->
        state + [output_metric: output.output]
      }
    )

    | yaml_to_json.run(
      key: "dataset_info",
      fromState: [ 
        "input": "input_dataset_info",
        "output": "output_dataset_info"
      ],
      toState: { id, output, state ->
        state + [output_dataset: output.output]
      }
    )

    | yaml_to_json.run(
      key: "task_info",
      fromState: [ 
        "input": "input_task_info",
        "output": "output_task_info"
      ],
      toState: { id, output, state ->
        state + [output_task: output.output]
      }
    )

    | get_results.run(
      fromState: [ 
        "input_scores": "input_scores",
        "input_execution" : "input_execution",
        "output": "output_scores"
      ],
      toState: { id, output, state ->
        state + [output_results: output.output]
      }
    )

    | map{ id, state ->
      def _meta = [join_id: id]

      def new_state = [
        output_scores: state.output_results,
        output_method_info: state.output_method,
        output_metric_info: state.output_metric,
        output_dataset_info: state.output_dataset,
        output_task_info: state.output_task,
        _meta: _meta
      ]

      ["output", new_state]
    }

  emit:
  output_ch
}