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

  // construct list of metrics
  metrics = [
    correlation,
    mse
  ]

  /****************************
   * EXTRACT DATASET METADATA *
   ****************************/
  dataset_ch = input_ch

    // store original id for later use
    | map{ id, state ->
      [id, state + ["_meta": [join_id: id]]]
    }

    // extract the dataset metadata
    | extract_metadata.run(
      key: "metadata_mod1",
      fromState: [input: "input_train_mod1"],
      toState: { id, output, state ->
        state + [
          dataset_uns_mod1: readYaml(output.output).uns
        ]
      }
    )

    | extract_metadata.run(
      key: "metadata_mod2",
      fromState: [input: "input_test_mod2"],
      toState: { id, output, state ->
        state + [
          dataset_uns_mod2: readYaml(output.output).uns
        ]
      }
    )

    | map{ id, state ->
      def rna_norm = state.dataset_uns_mod1.modality == "GEX" ? state.dataset_uns_mod1.normalization_id : state.dataset_uns_mod2.normalization_id
      [id, state + [rna_norm: rna_norm]]
    }

  /***************************
   * RUN METHODS AND METRICS *
   ***************************/
  score_ch = dataset_ch

    // run all methods
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        def norm = state.rna_norm
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
      fromState: { id, state, comp ->
        def new_args = [
          input_train_mod1: state.input_train_mod1,
          input_train_mod2: state.input_train_mod2,
          input_test_mod1: state.input_test_mod1
        ]
        if (comp.config.functionality.info.type == "control_method") {
          new_args.input_test_mod2 = state.input_test_mod2
        }
        new_args
      },

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
        input_test_mod2: "input_test_mod2", 
        input_prediction: "method_output"
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

  // extract the dataset metadata
  dataset_meta_ch = dataset_ch
    // only keep one of the normalization methods
    | filter{ id, state ->
      state.rna_norm == "log_cp10k"
    }
    | joinStates { ids, states ->
      // store the dataset metadata in a file
      def dataset_uns = states.collect{state ->
        def uns = state.dataset_uns_mod2.clone()
        uns.remove("normalization_id")
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