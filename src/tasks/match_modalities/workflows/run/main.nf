sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// import control methods
include { random_features }      from "$targetDir/match_modalities/control_methods/random_features/main.nf" params(params)
include { true_features }        from "$targetDir/match_modalities/control_methods/true_features/main.nf"   params(params)

// import methods
include { fastmnn }                  from "$targetDir/match_modalities/methods/fastmnn/main.nf"                 params(params)
include { scot }                 from "$targetDir/match_modalities/methods/scot/main.nf"                params(params)
include { harmonic_alignment }   from "$targetDir/match_modalities/methods/harmonic_alignment/main.nf"  params(params)
include { procrustes }           from "$targetDir/match_modalities/methods/procrustes/main.nf"          params(params)

// import metrics
include { knn_auc }              from "$targetDir/match_modalities/metrics/knn_auc/main.nf"             params(params)
include { mse }                  from "$targetDir/match_modalities/metrics/mse/main.nf"                 params(params)

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"                         params(params)

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { runComponents; joinStates; initializeTracer; writeJon; getPublishDir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

// read in pipeline config
config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

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

workflow {
    helpMessage(config)

    // create channel from input parameters with
    // arguments as defined in the config
    channelFromParams(params, config)
        | run_wf
}

    // run the workflow
workflow run_wf {
    take:
    input_ch

    main:
    output_ch = input_ch

    // based on the config file (config.vsh.yaml), run assertions on parameter sets
    // and fill in default values
    | preprocessInputs(config: config)

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
    | runComponents(
        components: methods,

        // // use the 'filter' argument to only run a method on the normalisation the component is asking for
        // filter: { id, state, config ->
        // def norm = state.normalization_id
        // def pref = config.functionality.info.preferred_normalization
        // // if the preferred normalisation is none at all,
        // // we can pass whichever dataset we want
        // (norm == "log_cp10k" && pref == "counts") || norm == pref
        // },

        // define a new 'id' by appending the method name to the dataset id
        id: { id, state, config ->
            id + "." + config.functionality.name
        },

        // use 'fromState' to fetch the arguments the component requires from the overall state
        fromState: { id, state, config ->
            def new_args = [
            input_mod1: state.input_mod1,
            input_mod2: state.input_mod2
            ]
            new_args
        },

        // use 'toState' to publish that component's outputs to the overall state
        toState: { id, output, state, config ->
            state + [
            method_id: config.functionality.name,
            method_output_mod1: output.output_mod1,
            method_output_mod2: output.output_mod2
            ]
        }
    )

        // run all metrics
    | runComponents(
        components: metrics,
        // use 'fromState' to fetch the arguments the component requires from the overall state
        fromState: [
            input_mod1: "method_output_mod1",
            input_mod2: "method_output_mod2"
        ],
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

    writeJsontraces, file("$publish_dir/traces.json"))
    // todo: add datasets logging
    writeJsonmethods.collect{it.config}, file("$publish_dir/methods.json"))
    writeJsonmetrics.collect{it.config}, file("$publish_dir/metrics.json"))
}
