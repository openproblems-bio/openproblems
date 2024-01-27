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

  // construct list of methods
  methods = [
    no_denoising,
    perfect_denoising,
    alra,
    dca,
    knn_smoothing,
    magic
  ]

  // construct list of metrics
  metrics = [
    mse,
    poisson
  ]

  /****************************
   * EXTRACT DATASET METADATA *
   ****************************/
  dataset_ch = input_ch
    // store join id
    | map{ id, state -> 
      [id, state + ["_meta": [join_id: id]]]
    }

    // extract the dataset metadata
    | extract_metadata.run(
      fromState: [input: "input_test"],
      toState: { id, output, state ->
        state + [
          dataset_uns: readYaml(output.output).uns
        ]
      }
    )
    
  /***************************
   * RUN METHODS AND METRICS *
   ***************************/
  score_ch = dataset_ch

    // run all methods
    | runEach(
      components: methods,
      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_train: "input_train",
        input_test: "input_test"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.functionality.name,
          method_output: output.output
        ]
      }
    )

    // run all metrics
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_test: "input_test", 
        input_denoised: "method_output"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.functionality.name,
          metric_output: output.output
        ]
      }
    )

  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  // TODO: can we store everything below in a separate helper function?
  // NOTE: the 'denoising' task doesn't use normalized data,
  // so code related to normalization_ids is commented out

  // extract the dataset metadata
  dataset_meta_ch = dataset_ch
    // // only keep one of the normalization methods
    // | filter{ id, state ->
    //   state.dataset_uns.normalization_id == "log_cp10k"
    // }
    | joinStates { ids, states ->
      // store the dataset metadata in a file
      def dataset_uns = states.collect{state ->
        def uns = state.dataset_uns.clone()
        // uns.remove("normalization_id")
        uns
      }
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      ["output", [output_dataset_info: dataset_uns_file]]
    }

  output_ch = score_ch

    // extract the scores
    | extract_metadata.run(
      key: "extract_scores",
      fromState: [input: "metric_output"],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      // store the method configs in a file
      def method_configs = methods.collect{it.config}
      def method_configs_yaml_blob = toYamlBlob(method_configs)
      def method_configs_file = tempFile("method_configs.yaml")
      method_configs_file.write(method_configs_yaml_blob)

      // store the metric configs in a file
      def metric_configs = metrics.collect{it.config}
      def metric_configs_yaml_blob = toYamlBlob(metric_configs)
      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(metric_configs_yaml_blob)

      def task_info_file = meta.resources_dir.resolve("task_info.yaml")

      // store the scores in a file
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      def new_state = [
        output_method_configs: method_configs_file,
        output_metric_configs: metric_configs_file,
        output_task_info: task_info_file,
        output_scores: score_uns_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

    // merge all of the output data 
    | mix(dataset_meta_ch)
    | joinStates{ ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}