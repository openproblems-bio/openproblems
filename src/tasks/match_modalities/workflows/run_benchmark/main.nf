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
    random_features,
    true_features,
    scot,
    harmonic_alignment,
    fastmnn,
    procrustes
  ]

  // collect metric list
  metrics = [
    knn_auc,
    mse
  ]

  output_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [ "input": "input_mod1" ],
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
        id + "." + comp.config.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: { id, state, comp ->
        def new_args = [
          input_mod1: state.input_mod1,
          input_mod2: state.input_mod2
        ]
        if (comp.config.functionality.info.type == "control_method") {
          new_args.input_solution_mod1 = state.input_solution_mod1
          new_args.input_solution_mod2 = state.input_solution_mod2
        }
        new_args
      },

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.functionality.name,
          method_output_mod1: output.output_mod1,
          method_output_mod2: output.output_mod2
        ]
      }
    )

      // run all metrics
    | runEach(
      components: metrics,
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_mod1: "method_output_mod1",
        input_mod2: "method_output_mod2",
        input_solution_mod1: "input_solution_mod1",
        input_solution_mod2: "input_solution_mod2"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.functionality.name,
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
