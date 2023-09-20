nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// import methods
include { bbknn } from "$targetDir/batch_integration/methods/bbknn/main.nf"
include { combat } from "$targetDir/batch_integration/methods/combat/main.nf"
include { scanorama_embed } from "$targetDir/batch_integration/methods/scanorama_embed/main.nf"
include { scanorama_feature } from "$targetDir/batch_integration/methods/scanorama_feature/main.nf"
include { scvi } from "$targetDir/batch_integration/methods/scvi/main.nf"

// import control methods
include { no_integration_batch } from "$targetDir/batch_integration/control_methods/no_integration_batch/main.nf"
include { random_embed_cell } from "$targetDir/batch_integration/control_methods/random_embed_cell/main.nf"
include { random_embed_cell_jitter } from "$targetDir/batch_integration/control_methods/random_embed_cell_jitter/main.nf"
include { random_integration } from "$targetDir/batch_integration/control_methods/random_integration/main.nf"

// import transformers
include { feature_to_embed } from "$targetDir/batch_integration/transformers/feature_to_embed/main.nf"
include { embed_to_graph } from "$targetDir/batch_integration/transformers/embed_to_graph/main.nf"

// import metrics
include { asw_batch } from "$targetDir/batch_integration/metrics/asw_batch/main.nf"
include { asw_label } from "$targetDir/batch_integration/metrics/asw_label/main.nf"
include { cell_cycle_conservation } from "$targetDir/batch_integration/metrics/cell_cycle_conservation/main.nf"
include { clustering_overlap } from "$targetDir/batch_integration/metrics/clustering_overlap/main.nf"
include { graph_connectivity } from "$targetDir/batch_integration/metrics/graph_connectivity/main.nf"
include { hvg_overlap } from "$targetDir/batch_integration/metrics/hvg_overlap/main.nf"
include { isolated_label_asw } from "$targetDir/batch_integration/metrics/isolated_label_asw/main.nf"
include { isolated_label_f1 } from "$targetDir/batch_integration/metrics/isolated_label_f1/main.nf"
include { kbet } from "$targetDir/batch_integration/metrics/kbet/main.nf"
include { lisi } from "$targetDir/batch_integration/metrics/lisi/main.nf"
include { pcr } from "$targetDir/batch_integration/metrics/pcr/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs; readYaml } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { publishState; runComponents; joinStates; initializeTracer; writeJson; getPublishDir; autoDetectStates } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

// collect method list
methods = [
  bbknn,
  combat,
  scanorama_embed,
  scanorama_feature,
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


workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
    | publishState([:])
}

workflow auto {
  autoDetectStates(params, config)
    | run_wf
    | publishState([:])
}

workflow run_wf {
  take:
  input_ch

  main:

  // process input parameter channel
  dataset_ch = input_ch
    | preprocessInputs(config: config)

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

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
  // todo: add datasets logging
  // writeJson(methods.collect{it.config}, file("$publish_dir/methods.json"))
  // writeJson(metrics.collect{it.config}, file("$publish_dir/metrics.json"))
}