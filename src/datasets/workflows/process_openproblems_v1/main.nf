nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// dataset loaders
include { openproblems_v1 } from "$targetDir/datasets/loaders/openproblems_v1/main.nf"

// normalization methods
include { log_cpm } from "$targetDir/datasets/normalization/log_cpm/main.nf"
include { log_scran_pooling } from "$targetDir/datasets/normalization/log_scran_pooling/main.nf"
include { sqrt_cpm } from "$targetDir/datasets/normalization/sqrt_cpm/main.nf"
include { l1_sqrt } from "$targetDir/datasets/normalization/l1_sqrt/main.nf"

// dataset processors
include { pca } from "$targetDir/datasets/processors/pca/main.nf"
include { hvg } from "$targetDir/datasets/processors/hvg/main.nf"
include { knn } from "$targetDir/datasets/processors/knn/main.nf"
include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { run_components; join_states; initialize_tracer; write_json; get_publish_dir } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initialize_tracer()

// normalization_methods = [log_cpm, log_scran_pooling, sqrt_cpm, l1_sqrt
normalization_methods = [log_cpm, sqrt_cpm, l1_sqrt]

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs(config: config)

    // fetch data from legacy openproblems
    | run_components(
      components: openproblems_v1,
      from_state: [
        "dataset_id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts", "sparse",
        "dataset_name", "data_url", "data_reference", "dataset_summary", "dataset_description", "dataset_organism"
      ],
      to_state: [ dataset: "output" ]
    )

    // run normalization methods
    | run_components(
      components: normalization_methods,
      id: { id, state, config -> id + "/" + config.functionality.name },
      from_state: [ input: "dataset" ],
      to_state: [
        normalization_id: config.functionality.name,
        output_normalization: "output"
      ]
    )

    | run_components(
      components: pca,
      from_state: [ input: "output_normalization" ],
      to_state: [ pca: "output" ]
    )

    | run_components(
      components: hvg,
      from_state: [ input: "pca" ],
      to_state: [ hvg: "output" ]
    )

    | run_components(
      components: knn,
      from_state: [ input: "hvg" ],
      to_state: [ knn: "output" ]
    )

    | run_components(
      components: check_dataset_schema,
      from_state: {id, state, config ->
        [
          input: state.knn,
          meta: state.output_meta,
          output: state.output_dataset,
          checks: null
        ]
      },
      to_state: [],
      auto: [publish: true]
    )

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = get_publish_dir()

  write_json(traces, file("$publish_dir/traces.json"))
  write_json(normalization_methods.collect{it.config}, file("$publish_dir/normalization_methods.json"))
}