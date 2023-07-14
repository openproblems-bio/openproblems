sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { mean_per_gene } from "$targetDir/predict_modality/control_methods/mean_per_gene/main.nf"
include { random_predict } from "$targetDir/predict_modality/control_methods/random_predict/main.nf"
include { zeros } from "$targetDir/predict_modality/control_methods/zeros/main.nf"
include { solution } from "$targetDir/predict_modality/control_methods/solution/main.nf"

// import methods
include { knnr_py } from "$targetDir/predict_modality/methods/knnr_py/main.nf"
include { knnr_r } from "$targetDir/predict_modality/methods/knnr_r/main.nf"
include { lm } from "$targetDir/predict_modality/methods/lm/main.nf"
include { newwave_knnr } from "$targetDir/predict_modality/methods/newwave_knnr/main.nf"
include { random_forest } from "$targetDir/predict_modality/methods/random_forest/main.nf"


// import metrics
include { correlation } from "$targetDir/predict_modality/metrics/correlation/main.nf"
include { mse } from "$targetDir/predict_modality/metrics/mse/main.nf"

// tsv generation component
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

    // run all methods
    | run_components(
      components: methods,

      // // use the 'filter' argument to only run a method on the normalisation the component is asking for
      // filter: { id, state, config ->
      //   def norm = state.normalization_id
      //   def pref = config.functionality.info.preferred_normalization
      //   // if the preferred normalisation is none at all,
      //   // we can pass whichever dataset we want
      //   (norm == "log_cpm" && pref == "counts") || norm == pref
      // },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, config ->
        id + "." + config.functionality.name
      },

      // use 'from_state' to fetch the arguments the component requires from the overall state
      from_state: { id, state, config ->
        def new_args = [
          input_train_mod1: state.input_train_mod1,
          input_train_mod2: state.input_train_mod2,
          input_test_mod1: state.input_test_mod1
        ]
        if (config.functionality.info.type == "control_method") {
          new_args.input_test_mod2 = state.input_test_mod2
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
      from_state: [
        input_test_mod2: "input_test_mod2", 
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