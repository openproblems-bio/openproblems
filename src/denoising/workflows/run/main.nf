nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { no_denoising } from "$targetDir/denoising/control_methods/no_denoising/main.nf"
include { perfect_denoising } from "$targetDir/denoising/control_methods/perfect_denoising/main.nf"

// import methods
include { magic } from "$targetDir/denoising/methods/magic/main.nf"


// import metrics
include { mse } from "$targetDir/denoising/metrics/mse/main.nf"
include { poisson } from "$targetDir/denoising/metrics/poisson/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
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
      method: ["input_train", "input_test"],
      metric: ["input_test"],
      output: ["output"]
    )

    // run methods
    | getWorkflowArguments(key: "method")
    | (
      no_denoising &
      perfect_denoising & 
      magic
    )
    | mix

    // construct tuples for metrics
    | pmap{ id, file, passthrough ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      // combine prediction with solution
      def newData = [ input_denoised: file, input_test: passthrough.metric.input_test ]
      [ newId, newData, passthrough ]
    }
    
    // run metrics
    | (mse & poisson)
    | mix

    // convert to tsv  
    | toSortedList
    | map{ it -> [ "combined", it.collect{ it[1] } ] + it[0].drop(2) }
    | getWorkflowArguments(key: "output")
    | extract_scores.run(
        auto: [ publish: true ]
    )

  emit:
  output_ch
}