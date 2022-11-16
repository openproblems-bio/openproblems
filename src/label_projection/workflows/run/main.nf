nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { true_labels } from "$targetDir/label_projection/control_methods/true_labels/main.nf"
include { majority_vote } from "$targetDir/label_projection/control_methods/majority_vote/main.nf"
include { random_labels } from "$targetDir/label_projection/control_methods/random_labels/main.nf"

// import methods
include { knn_classifier } from "$targetDir/label_projection/methods/knn_classifier/main.nf"
include { mlp } from "$targetDir/label_projection/methods/mlp/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
// include { scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scanvi_hvg/main.nf"
// include { scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scanvi_all_genes/main.nf"
// include { scarches_scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_all_genes/main.nf"
// include { scarches_scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_hvg/main.nf"

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
  output_ch = input_ch
    
    // split params for downstream components
    | setWorkflowArguments(
      method: ["input_train", "input_test"],
      metric: [ "input_solution" ]
    )

    // run method
    | getWorkflowArguments(key: "method")
    | majority_vote

    // run metric
    | getWorkflowArguments(key: "metric", inputKey: "input_prediction")
    | (accuracy & f1)
    | mix

  emit:
  output_ch
}
// workflow {
//     load_data
//         | randomize
//         | subsample.run(
//             map: { [it[0], [input: it[1], even: true]] }
//         )
//         | (log_cpm & log_scran_pooling)
//         | mix
//         | map { unique_file_name(it) }
//         | (knn_classifier & mlp0 & lr0 & random_labels & majority_vote & true_labels)
//         | mix
//         | map { unique_file_name(it) }
//         | (accuracy & f1a)
//         | mix
//         | toSortedList
//         | map{ it -> [ "combined", [ input: it.collect{ it[1] } ] ] }
//         | extract_scores.run(
//             args: [column_names: "dataset_id:normalization_method:method_id:metric_id:metric_value"],
//             auto: [ publish: true ]
//         )
// }
