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

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [ "input": "input_train_mod1" ],
      toState: { id, output, state ->
        state + (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
      }
    )

    | check_dataset_schema.run(
      fromState: [ "input": "input_train_mod2" ],
      toState: { id, output, state ->
        state + (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
      }
    )

    // run all methods
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        def norm = state.normalization_id
        def pref = comp.config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        (norm == "log_cp10k" && pref == "counts") || norm == pref
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
          input_train_mod1: "input_train_mod1",
          input_train_mod2: "input_train_mod2",
          input_test_mod1: "input_test_mod1",
          input_test_mod2: "input_test_mod2"
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
        input_test_mod2: "input_test_mod2", 
        input_prediction: "method_output"
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