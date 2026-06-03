include { runEach; joinStates; readYaml; toYamlBlob } from "${moduleDir}/workflowHelper.nf"

def run_methods(args) {
  // required args
  def methodComponents = args.methods
  def methodFromState = args.fromState
  def methodToState = args.toState

  assert methodComponents, "argument methods must be defined"
  assert methodFromState, "argument fromState must be defined"
  assert methodToState, "argument toState must be defined"

  // optional args
  def keyPrefix = args.keyPrefix ?: ""
  def methodFilter = args.filter
  def methodAuto = args.auto ?: [:]
  if (!methodFilter) {
    methodFilter = { id, state, comp ->
      !state.method_ids || state.method_ids.contains(comp.config.name)
    }
  }

  // add the key prefix to the method names
  if (keyPrefix && keyPrefix != "") {
    methodComponents = methodComponents.collect{ method ->
      method.run(key: keyPrefix + method.config.name)
    }
  }

  workflow methodsWf {
    take: input_ch

    main:
    output_ch = input_ch
      // run all methods
      | runEach(
        components: methodComponents,
        filter: methodFilter,
        id: { id, state, comp ->
          id + "." + comp.config.name
        },
        fromState: methodFromState,
        toState: methodToState,
        auto: methodAuto
      )

    emit: output_ch
  }
  return methodsWf
}

def run_metrics(args) {
  // required args
  def metricComponents = args.metrics
  def metricFromState = args.fromState
  def metricToState = args.toState

  assert metricComponents, "argument metrics must be defined"
  assert metricFromState, "argument fromState must be defined"
  assert metricToState, "argument toState must be defined"

  // optional args
  def keyPrefix = args.keyPrefix ?: ""
  def metricFilter = args.filter
  def metricAuto = args.auto ?: [:]
  if (!metricFilter) {
    metricFilter = { id, state, comp ->
      !state.metric_ids || state.metric_ids.contains(comp.config.name)
    }
  }

  // add the key prefix to the metric names
  if (keyPrefix && keyPrefix != "") {
    metricComponents = metricComponents.collect{ metric ->
      metric.run(key: keyPrefix + metric.config.name)
    }
  }

  workflow metricsWf {
    take: input_ch

    main:
    output_ch = input_ch
      // run all metrics
      | runEach(
        components: metricComponents,
        filter: metricFilter,
        id: { id, state, comp ->
          id + "." + comp.config.name
        },
        fromState: metricFromState,
        toState: metricToState,
        auto: metricAuto
      )

    emit: output_ch
  }
  return metricsWf
}

def extract_scores(args) {
  // required args
  def extract_uns_metadata_component = args.extract_uns_metadata_component

  assert extract_uns_metadata_component, "argument extract_uns_metadata_component must be defined"

  // optional args
  def keyPrefix = args.keyPrefix ?: ""

  workflow extractScoresWf {
    take: input_ch

    main:
    output_ch = input_ch
      | extract_uns_metadata_component.run(
        key: "${keyPrefix}score_uns",
        fromState: [input: "metric_output"],
        toState: { id, output, state ->
          def outputYaml = readYaml(output.output)
          if (!outputYaml.uns) {
            throw new Exception("id '$id': No uns found in provided metric output")
          }
          state + [score_uns: outputYaml.uns]
        }
      )
    | joinStates { ids, states ->
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)
      
      ["output", [output_scores: score_uns_file]]
    }

    emit: output_ch
  }
  return extractScoresWf
}


def create_metadata_files(args) {
  // required args
  def meta_ = args.meta
  def datasetFile_ = args.datasetFile
  def methods_ = args.methods
  def metrics_ = args.metrics
  def extract_uns_metadata_component = args.extract_uns_metadata_component

  assert meta_, "argument meta must be defined"
  assert methods_, "argument methods must be defined"
  assert metrics_, "argument metrics must be defined"
  assert datasetFile_, "argument datasetFile must be defined"
  assert extract_uns_metadata_component, "argument extract_uns_metadata_component must be defined"

  // optional args
  def datasetUnsModifier = args.datasetUnsModifier
  def datasetFilter = args.filter

  workflow createMetadataFilesWf {
    take: input_ch

    main:
    output_ch = input_ch

      | map{ id, state ->
        [id, state + ["_meta": [join_id: id]]]
      }

      | extract_uns_metadata_component.run(
        key: "dataset_uns",
        fromState: [input: datasetFile_],
        toState: { id, output, state ->
          def outputYaml = readYaml(output.output)
          if (!outputYaml.uns) {
            throw new Exception("id '$id': No uns found in provided dataset")
          }
          state + [ dataset_uns: outputYaml.uns ]
        }
      )

      | filter{ id, state ->
        if (datasetFilter) {
          datasetFilter(id, state)
        } else {
          true
        }
      }
    
      | joinStates { ids, states ->
        assert states.size() > 0, "no states found"
        assert states[0]._meta, "no _meta found in state[0]"
        assert states.every{it.dataset_uns}, "not all states have dataset_info"

        // combine the dataset info into one file
        def dataset_uns = states.collect{state ->
          def uns = state.dataset_uns
          if (datasetUnsModifier) {
            uns = datasetUnsModifier(uns.clone())
          }
          uns
        }
        def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
        def dataset_uns_file = tempFile("dataset_uns.yaml")
        dataset_uns_file.write(dataset_uns_yaml_blob)

        // store the method configs in a file
        def method_configs = methods_.collect{it.config}
        def method_configs_yaml_blob = toYamlBlob(method_configs)
        def method_configs_file = tempFile("method_configs.yaml")
        method_configs_file.write(method_configs_yaml_blob)

        // store the metric configs in a file
        def metric_configs = metrics_.collect{it.config}
        def metric_configs_yaml_blob = toYamlBlob(metric_configs)
        def metric_configs_file = tempFile("metric_configs.yaml")
        metric_configs_file.write(metric_configs_yaml_blob)

        def task_info_file = meta_.resources_dir.resolve("_viash.yaml")

        def new_state = [
          output_dataset_info: dataset_uns_file,
          output_method_configs: method_configs_file,
          output_metric_configs: metric_configs_file,
          output_task_info: task_info_file,
          _meta: states[0]._meta
        ]
        ["output", new_state]
      }
    emit: output_ch
  }
  return createMetadataFilesWf
}
