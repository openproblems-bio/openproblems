sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

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
include { rmse } from "$targetDir/dimensionality_reduction/metrics/rmse/main.nf"
include { trustworthiness } from "$targetDir/dimensionality_reduction/metrics/trustworthiness/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { run_components; join_states; initialize_tracer; write_json; get_publish_dir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initialize_tracer()

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
  rmse,
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

    // run all methods
    | run_components(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, config ->
        def norm = state.normalization_id
        def pref = config.functionality.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        (norm == "log_cpm" && pref == "counts") || norm == pref
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, config ->
        id + "." + config.functionality.name
      },

      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: { id, state, config ->
        def new_args = [
          input: state.input
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_solution = state.input_solution
        }
        new_args
      },

      // use 'to_state' to publish that component's outputs to the overall state
      to_state: { id, output, config ->
        [
          method_id: config.functionality.name,
          method_output: output.output
        ]
      }
    )

    // run all metrics
    | run_components(
      components: metrics,
      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: { id, state, config ->
        [
          input_solution: state.input_solution,
          input_embedding: state.method_output
        ]
      },
      // use 'to_state' to publish that component's outputs to the overall state
      to_state: { id, output, config ->
        [
          metric_id: config.functionality.name,
          metric_output: output.output
        ]
      }
    )

    // join all events into a new event where the new id is simply "output" and the new state consists of:
    //   - "input": a list of score h5ads
    //   - "output": the output argument of this workflow
    | join_states(
      apply: { ids, states ->
        ["output", [
          input: states.collect{it.metric_output},
          output: states[0].output
        ]]
      }
    )

    // convert to tsv and publish
    | extract_scores.run(
      auto: [publish: true]
    )

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = get_publish_dir()

  write_json(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  write_json(methods.collect{it.config}, file("$publish_dir/methods.json"))
  write_json(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}