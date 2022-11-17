nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

include { openproblems_v1 } from "$targetDir/datasets/loaders/openproblems_v1/main.nf"
include { log_cpm } from "$targetDir/datasets/normalization/log_cpm/main.nf"
include { log_scran_pooling } from "$targetDir/datasets/normalization/log_scran_pooling/main.nf"
include { sqrt_cpm } from "$targetDir/datasets/normalization/sqrt_cpm/main.nf"

include { readConfig; viashChannel; helpMessage } from sourceDir + "/nxf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/nxf_utils/DataFlowHelper.nf"

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
      loader: ["id", "obs_celltype", "obs_batch", "obs_tissue", "layer_counts", "sparse"],
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

    // run scran normalisation
    | log_scran_pooling

    // run sqrt normalisation and publish
    | getWorkflowArguments(key: "sqrt_cpm")
    | sqrt_cpm.run(
      auto: [ publish: true ]
    )

    // clean up channel
    | pmap{id, data, passthrough -> [id, data]}

  emit:
  output_ch
}