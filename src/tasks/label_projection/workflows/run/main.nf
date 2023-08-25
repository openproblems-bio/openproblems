sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { true_labels } from "$targetDir/label_projection/control_methods/true_labels/main.nf"
include { majority_vote } from "$targetDir/label_projection/control_methods/majority_vote/main.nf"
include { random_labels } from "$targetDir/label_projection/control_methods/random_labels/main.nf"

// import methods
include { knn } from "$targetDir/label_projection/methods/knn/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
include { mlp } from "$targetDir/label_projection/methods/mlp/main.nf"
include { scanvi } from "$targetDir/label_projection/methods/scanvi/main.nf"
include { scanvi_scarches } from "$targetDir/label_projection/methods/scanvi_scarches/main.nf"
// include { seurat_transferdata } from "$targetDir/label_projection/methods/seurat_transferdata/main.nf"
include { xgboost } from "$targetDir/label_projection/methods/xgboost/main.nf"

// import metrics
include { accuracy } from "$targetDir/label_projection/metrics/accuracy/main.nf"
include { f1 } from "$targetDir/label_projection/metrics/f1/main.nf"

// convert scores to tsv
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { run_components; join_states; initialize_tracer; write_json; get_publish_dir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initialize_tracer()

// collect method list
methods = [
  true_labels,
  majority_vote,
  random_labels,
  knn,
  logistic_regression,
  mlp,
  scanvi,
  scanvi_scarches,
  // seurat_transferdata,
  xgboost
]

// collect metric list
metrics = [
  accuracy,
  f1
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
    // based on the config file (config.vsh.yaml), run assertions on parameter sets
    // and fill in default values
    | preprocessInputs(config: config)

    | view{ id, state ->
      "input event: [id: $id, dataset_id: $state.dataset_id, normalization_id: $state.normalization_id, ...]"
    }

    // run all methods
    | run_components(
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

      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: { id, state, config ->
        def new_args = [
          input_train: state.input_train,
          input_test: state.input_test
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

    | view{ id, state ->
      "after method: [id: $id, dataset_id: $state.dataset_id, normalization_id: $state.normalization_id, method_id: $state.method_id, ...]"
    }

    // run all metrics
    | run_components(
      components: metrics,
      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: [
        input_solution: "input_solution", 
        input_prediction: "method_output"
      ],
      // use 'to_state' to publish that component's outputs to the overall state
      to_state: { id, output, config ->
        [
          metric_id: config.functionality.name,
          metric_output: output.output
        ]
      }
    )

    | view{ id, state ->
      "after metric: [id: $id, dataset_id: $state.dataset_id, normalization_id: $state.normalization_id, method_id: $state.method_id, metric_id: $state.metric_id, ...]"
    }

    // join all events into a new event where the new id is simply "output" and the new state consists of:
    //   - "input": a list of score h5ads
    //   - "output": the output argument of this workflow
    | join_states{ ids, states ->
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
  def publish_dir = get_publish_dir()

  write_json(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  write_json(methods.collect{it.config}, file("$publish_dir/methods.json"))
  write_json(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}