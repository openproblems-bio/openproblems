// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

workflow run_wf {
  take:
  input_ch

  main:

  // collect method list
  methods = [
    mean_per_gene,
    random_predict,
    zeros,
    solution,
    knnr_py,
    knnr_r,
    lm,
    newwave_knnr,
    random_forest
  ]

  // collect metric list
  metrics = [
    correlation,
    mse
  ]

  output_ch = input_ch

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [ "input": "input_train_mod1" ],
      toState: { id, output, state ->
        // load output yaml file
        def metadata = (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
        // add metadata from file to state
        state + metadata
      }
    )

    // run all methods
    | runComponents(
      components: methods,

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, config ->
        id + "." + config.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: { id, state, config ->
        def new_args = [
          input_train_mod1: state.input_train_mod1,
          input_train_mod2: state.input_train_mod2,
          input_test_mod1: state.input_test_mod1
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_test_mod2 = state.input_test_mod2
        }
        new_args
      },

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, config ->
        state + [
          method_id: config.functionality.name,
          method_output: output.output
        ]
      }
    )

    // run all metrics
    | runComponents(
      components: metrics,
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_test_mod2: "input_test_mod2", 
        input_prediction: "method_output"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, config ->
        state + [
          metric_id: config.functionality.name,
          metric_output: output.output
        ]
      }
    )

    // join all events into a new event where the new id is simply "output" and the new state consists of:
    //   - "input": a list of score h5ads
    //   - "output": the output argument of this workflow
    | joinStates{ ids, states ->
      def new_id = "output"
      def new_state = [
        input: states.collect{it.metric_output},
        output: states[0].output
      ]
      [new_id, new_state]
    }

    // convert to tsv and publish
    | extract_scores.run(
      auto: [publish: true]
    )

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  writeJson(methods.collect{it.config}, file("$publish_dir/methods.json"))
  writeJson(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}