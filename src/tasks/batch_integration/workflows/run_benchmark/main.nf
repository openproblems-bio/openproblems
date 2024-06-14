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
    no_integration_batch_embed,
    no_integration_global_embed,
    no_integration_global_feature,
    no_integration_global_graph,
    perfect_integration_celltype_embed,
    perfect_integration_celltype_jitter_embed,
    random_integration_batch_embed,
    random_integration_batch_feature,
    random_integration_batch_graph,
    random_integration_celltype_embed,
    random_integration_celltype_feature,
    random_integration_celltype_graph,
    random_integration_global_embed,
    random_integration_global_feature,
    random_integration_global_graph,
  ]

  // construct list of metrics
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
    pcr
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
      fromState: [input: "input_solution"],
      toState: { id, output, state ->
        state + [
          dataset_uns: readYaml(output.output).uns
        ]
      }
    )

  /***************************
   * RUN METHODS AND METRICS *
   ***************************/
  // run all methods
  method_out_ch1 = dataset_ch
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        def norm = state.dataset_uns.normalization_id
        def pref = comp.config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        def norm_check = (norm == "log_cp10k" && pref == "counts") || norm == pref
        def method_check = !state.method_ids || state.method_ids.contains(comp.config.functionality.name)

        method_check && norm_check
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
      id: { id, state, comp ->
        id + "_f2e"
      },
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
      id: { id, state, comp ->
        id + "_e2g"
      },
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
  score_ch = method_out_ch3
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },
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


  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  // TODO: can we store everything below in a separate helper function?

  // extract the dataset metadata
  dataset_meta_ch = dataset_ch
    // only keep one of the normalization methods
    | filter{ id, state ->
      state.dataset_uns.normalization_id == "log_cp10k"
    }
    | joinStates { ids, states ->
      // store the dataset metadata in a file
      def dataset_uns = states.collect{state ->
        def uns = state.dataset_uns.clone()
        uns.remove("normalization_id")
        uns
      }
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      ["output", [output_dataset_info: dataset_uns_file]]
    }

  output_ch = score_ch
    // extract scores
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

      // store the task info in a file
      def task_info_file = meta.resources_dir.resolve("task_info.yaml")

      // store the scores in a file
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      // create state
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
