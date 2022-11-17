nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "target/nextflow"

// import methods
include { baseline } from "$targetDir/denoising/methods/baseline/main_nf"
include { magic } from "$targetDir/denoising/methods/magic/main_nf"


// import metrics
include { mse } from "$targetDir/denoising/metrics/mse/main_nf"
include { poisson } from "$targetDir/denoising/metrics/poisson/main_nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/nxf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/nxf_utils/DataFlowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}