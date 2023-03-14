nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { openproblems_v1 } from "$targetDir/datasets/loaders/openproblems_v1/main.nf"
include { log_cpm } from "$targetDir/datasets/normalization/log_cpm/main.nf"
include { log_scran_pooling } from "$targetDir/datasets/normalization/log_scran_pooling/main.nf"
include { sqrt_cpm } from "$targetDir/datasets/normalization/sqrt_cpm/main.nf"
include { pca } from "$targetDir/datasets/processors/pca/main.nf"
include { hvg } from "$targetDir/datasets/processors/hvg/main.nf"
include { knn } from "$targetDir/datasets/processors/knn/main.nf"

include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // split params for downstream components
    | setWorkflowArguments(
      loader: [
        "dataset_id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts", "sparse",
        "dataset_name", "data_url", "data_reference", "dataset_summary", "dataset_description", "dataset_organism"
      ],
      output: [ "output" ]
    )

    // fetch data from legacy openproblems
    | getWorkflowArguments(key: "loader")
    | openproblems_v1

    // run normalization methods
    | (log_cpm & log_scran_pooling & sqrt_cpm)
    | mix

    // make id unique again
    | pmap{ id, file ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      [ newId, file ]
    }

    | pca
    | hvg

    | getWorkflowArguments(key: "output")
    | knn.run(
      auto: [ publish: true ]
    )

    // clean up channel
    | pmap{id, data, passthrough -> [id, data]}

  emit:
  output_ch
}