nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { true_labels } from "$targetDir/label_projection/control_methods/true_labels/main.nf"
include { majority_vote } from "$targetDir/label_projection/control_methods/majority_vote/main.nf"
include { random_labels } from "$targetDir/label_projection/control_methods/random_labels/main.nf"

// import methods
include { knn } from "$targetDir/label_projection/methods/knn/main.nf"
include { mlp } from "$targetDir/label_projection/methods/mlp/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
include { scanvi } from "$targetDir/label_projection/methods/scanvi/main.nf"
include { seurat_transferdata } from "$targetDir/label_projection/methods/seurat_transferdata/main.nf"
include { xgboost } from "$targetDir/label_projection/methods/xgboost/main.nf"

// import metrics
include { accuracy } from "$targetDir/label_projection/metrics/accuracy/main.nf"
include { f1 } from "$targetDir/label_projection/metrics/f1/main.nf"

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

/** 
  * Helper function for making sure the right data gets passed to a method
  */
def configureMethod = { module ->
  // if method is a control method, add the solution to the data field
  def methodType = module.config.functionality.info.type
  def mapFun = null
  if (methodType != "method") {
    mapFun = { tup ->
      out = tup.clone()
      out[1] = out[1] + [input_solution: out[2].metric.input_solution]
      out
    }
  }

  // only let datasets with the right preferred normalization through
  // if preferred is 'counts', use any normalization method
  def preferred = module.config.functionality.info.preferred_normalization
  if (preferred == "counts") {
    preferred = "log_cpm"
  }

  def filterFun = { tup ->
    tup[1].normalization_id == preferred
  }

  // return module
  module.run(
    map: mapFun,
    filter: filterFun
  )
}

workflow run_wf {
  take:
  input_ch

  main:

  output_ch = input_ch
    
    // split params for downstream components
    | setWorkflowArguments(
      method: ["input_train", "input_test", "normalization_id", "dataset_id"],
      metric: ["input_solution"],
      output: ["output"]
    )

    // run methods
    | getWorkflowArguments(key: "method")
    | (
      configureMethod(true_labels) &
      configureMethod(random_labels) &
      configureMethod(majority_vote) &
      configureMethod(knn) &
      configureMethod(logistic_regression) &
      configureMethod(mlp) &
      configureMethod(scanvi) &
      configureMethod(seurat_transferdata) &
      configureMethod(xgboost)
    )
    | mix

    // construct tuples for metrics
    | pmap{ id, file, passthrough ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      // combine prediction with solution
      def newData = [ input_prediction: file, input_solution: passthrough.metric.input_solution ]
      [ newId, newData, passthrough ]
    }

    // run metrics
    | (accuracy & f1)
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