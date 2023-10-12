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
  // construct a map of methods (id -> method_module)
  methods = [
    no_denoising,
    perfect_denoising,
    alra,
    dca,
    knn_smoothing,
    magic
  ]

  metrics = [
    mse,
    poisson
  ]
  

  output_ch = input_ch
    
    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [ "input": "input_train" ],
      toState: { id, output, state ->
        // load output yaml file
        def metadata = (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
        // add metadata from file to state
        state + metadata
      }
    )

    // run all methods
    | runEach(
      components: methods,

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
          input_train: "input_train",
          input_test: "input_test"
        ],


      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.functionality.name,
          method_output: output.output
        ]
      }
    )

    // run all metrics
    | runEach(
      components: metrics,
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_test: "input_test", 
        input_denoised: "method_output"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.functionality.name,
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
        _meta: states[0]._meta
      ]
      [new_id, new_state]
    }
    
    // convert to tsv and publish
  | extract_scores.run(
    fromState: ["input"],
    toState: ["output"]
  )

  | setState(["output", "_meta"])

  emit:
  output_ch
}