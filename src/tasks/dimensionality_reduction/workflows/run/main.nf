sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// import control methods
include { random_features } from "$targetDir/dimensionality_reduction/control_methods/random_features/main.nf"
include { true_features } from "$targetDir/dimensionality_reduction/control_methods/true_features/main.nf"

// import methods
include { densmap } from "$targetDir/dimensionality_reduction/methods/densmap/main.nf"
// include { ivis } from "$targetDir/dimensionality_reduction/methods/ivis/main.nf"
include { neuralee } from "$targetDir/dimensionality_reduction/methods/neuralee/main.nf"
include { pca } from "$targetDir/dimensionality_reduction/methods/pca/main.nf"
include { phate } from "$targetDir/dimensionality_reduction/methods/phate/main.nf"
include { tsne } from "$targetDir/dimensionality_reduction/methods/tsne/main.nf"
include { umap } from "$targetDir/dimensionality_reduction/methods/umap/main.nf"

// import metrics
include { coranking } from "$targetDir/dimensionality_reduction/metrics/coranking/main.nf"
include { density_preservation } from "$targetDir/dimensionality_reduction/metrics/density_preservation/main.nf"
include { distance_correlation } from "$targetDir/dimensionality_reduction/metrics/distance_correlation/main.nf"
include { trustworthiness } from "$targetDir/dimensionality_reduction/metrics/trustworthiness/main.nf"

// convert scores to tsv
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { runComponents; joinStates; initializeTracer; writeJson; getPublishDir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

// collect method list
methods = [
  random_features,
  true_features,
  densmap,
  neuralee,
  pca,
  phate,
  tsne,
  umap
]

// collect metric list
metrics = [
  coranking,
  density_preservation,
  distance_correlation,
  trustworthiness
]

workflow {
  helpMessage(config)

  // create channel from input parameters with
  // arguments as defined in the config
  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs(config: config)

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [input: "input_dataset"],
      toState: { id, output, state ->
        state + (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
      }
    )

    // run all methods
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
      fromState: { id, state, config ->
        def new_args = [
          input: state.input_dataset
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_solution = state.input_solution
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
      fromState: { id, state, config ->
        [
          input_solution: state.input_solution,
          input_embedding: state.method_output
        ]
      },
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