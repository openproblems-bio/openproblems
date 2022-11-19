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
// include { scarches } from "$targetDir/label_projection/methods/scarches/main.nf"

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


workflow run_wf {
  take:
  input_ch

  main:
  def addSolution = { tup ->
    out = tup.clone()
    out[1] = out[1] + [input_solution: out[2].metric.input_solution]
    out
  }

  output_ch = input_ch
    | filter{it[1].normalization_id == "log_cpm"}
    
    // split params for downstream components
    | setWorkflowArguments(
      method: ["input_train", "input_test", "normalization_id", "dataset_id"],
      metric: ["input_solution"],
      output: ["output"]
    )

    // run methods
    // TODO: these filters don't work atm.
    | getWorkflowArguments(key: "method")
    | (
      true_labels.run(map: addSolution, filter: {it[1].normalization_id == "log_cpm"}) & 
      random_labels.run(filter: {it[1].normalization_id == "log_cpm"}) & 
      majority_vote.run(filter: {it[1].normalization_id == "log_cpm"}) & 
      knn.run(filter: {it[1].normalization_id == "log_cpm"}) & 
      logistic_regression.run(filter: {it[1].normalization_id == "log_cpm"}) &
      mlp.run(filter: {it[1].normalization_id == "log_cpm"}) &
      scanvi.run(filter: {it[1].normalization_id == "log_cpm"})
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