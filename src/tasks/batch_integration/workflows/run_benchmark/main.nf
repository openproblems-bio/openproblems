// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = collectTraces()

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

    | map { id, state -> 
      def newId = id.replaceAll(/\//, "_")

      [newId, state]
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
    | runComponents(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, config ->
        def norm = state.normalization_id
        def pref = config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        (norm == "log_cp10k" && pref == "counts") || norm == pref
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, config ->
        id + "." + config.functionality.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [input: "input_dataset"],

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, config ->
        state + [
          method_id: config.functionality.name,
          method_output: output.output,
          method_subtype: config.functionality.info.subtype
        ]
      }
    )
  
  // append feature->embed transformations
  method_out_ch2 = method_out_ch1
    | runComponents(
      components: feature_to_embed,
      filter: { id, state, config -> state.method_subtype == "feature"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, config ->
        state + [
          method_output: output.output,
          method_subtype: config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch1)

  // append embed->graph transformations
  method_out_ch3 = method_out_ch2
    | runComponents(
      components: embed_to_graph,
      filter: { id, state, config -> state.method_subtype == "embedding"},
      fromState: [ input: "method_output" ],
      toState: { id, output, state, config ->
        state + [
          method_output: output.output,
          method_subtype: config.functionality.info.subtype
        ]
      }
    )
    | mix(method_out_ch2)

  // run metrics
  output_ch = method_out_ch3
    | runComponents(
      components: metrics,
      filter: { id, state, config ->
        state.method_subtype == config.functionality.info.subtype
      },
      fromState: [
        input_integrated: "method_output",
        input_solution: "input_solution"
      ],
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
  | extract_scores

  emit:
  output_ch
}

workflow auto {
  findStates(params, thisConfig)
    | run_wf
    | publishStates([key: thisConfig.functionality.name])
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  // writeJson(methods.collect{it.config}, file("$publish_dir/methods.json"))
  // writeJson(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}

def joinStates(Closure apply_) {
  workflow joinStatesWf {
    take: input_ch
    main:
    output_ch = input_ch
      | toSortedList
      | filter{ it.size() > 0 }
      | map{ tups ->
        def ids = tups.collect{it[0]}
        def states = tups.collect{it[1]}
        apply_(ids, states)
      }

    emit: output_ch
  }
  return joinStatesWf
}