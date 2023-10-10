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
    bbknn,
    combat,
    fastmnn_embedding,
    fastmnn_feature,
    liger,
    mnn_correct,
    mnnpy,
    pyliger,
    scalex_embed,
    scalex_feature,
    scanorama_embed,
    scanorama_feature,
    scanvi,
    scvi,
    no_integration_batch,
    random_embed_cell,
    random_embed_cell_jitter,
    random_integration
  ]

  // collect metric list
  metrics = [
    asw_batch,
    asw_label,
    cell_cycle_conservation,
    clustering_overlap,
    graph_connectivity,
    hvg_overlap,
    isolated_label_asw,
    isolated_label_f1,
    kbet,
    lisi,
    pcr
  ]

  // process input parameter channel
  dataset_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [input: "input_dataset"],
      toState: { id, output, state ->
        state + (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
      }
    )

    // run all methods
  method_out_ch1 = dataset_ch
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
      fromState: [input: "input_dataset"],

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.functionality.name,
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
  

  // append feature->embed transformations
  method_out_ch2 = method_out_ch1
    | runEach(
      components: feature_to_embed,
      filter: { id, state, comp -> state.method_subtype == "feature"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, comp ->
        state + [
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch1)

  // append embed->graph transformations
  method_out_ch3 = method_out_ch2
    | runEach(
      components: embed_to_graph,
      filter: { id, state, comp -> state.method_subtype == "embedding"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, comp ->
        state + [
          method_output: output.output,
          method_subtype: comp.config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch2)

  // run metrics
  output_ch = method_out_ch3
    | runEach(
      components: metrics,
      filter: { id, state, comp ->
        state.method_subtype == comp.config.functionality.info.subtype
      },
      fromState: [
        input_integrated: "method_output",
        input_solution: "input_solution"
      ],
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
