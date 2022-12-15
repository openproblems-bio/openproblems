nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { high_dim_pca } from "$targetDir/dimensionality_reduction/control_methods/high_dim_pca/main.nf"
include { random_features } from "$targetDir/dimensionality_reduction/control_methods/random_features/main.nf"

// import methods
include { umap } from "$targetDir/dimensionality_reduction/methods/umap/main.nf"
include { densmap } from "$targetDir/dimensionality_reduction/methods/densmap/main.nf"
include { phate } from "$targetDir/dimensionality_reduction/methods/phate/main.nf"
include { tsne } from "$targetDir/dimensionality_reduction/methods/tsne/main.nf"
include { pca } from "$targetDir/dimensionality_reduction/methods/pca/main.nf"
include { neuralee } from "$targetDir/dimensionality_reduction/methods/neuralee/main.nf"

// import metrics
include { rmse } from "$targetDir/dimensionality_reduction/metrics/rmse/main.nf"
include { trustworthiness } from "$targetDir/dimensionality_reduction/metrics/trustworthiness/main.nf"
include { density } from "$targetDir/dimensionality_reduction/metrics/density/main.nf"
include { nn_ranking } from "$targetDir/dimensionality_reduction/metrics/nn_ranking/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFilter as pfilter } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// construct a map of methods (id -> method_module)
methods = [ random_features, high_dim_pca, umap, densmap, phate, tsne, pca, neuralee ]
  .collectEntries{method ->
    [method.config.functionality.name, method]
  }

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
    | view{"step 0: $it"}
    | setWorkflowArguments(
      preprocess: ["dataset_id", "normalization_id"],
      method: ["input"],
      metric: ["input_test"],
      output: ["output"]
    )
    // multiply events by the number of method
    | getWorkflowArguments(key: "preprocess")
    | add_methods

    // filter the normalization methods that a method actually prefers
    | check_filtered_normalization_id
    // add input_solution to data for the positive controls
    | controls_can_cheat

    // run methods
    | getWorkflowArguments(key: "method")
    | run_methods

    // run metrics
    | getWorkflowArguments(key: "metric", inputKey: "input_reduced")
    | run_metrics

    // convert to tsv  
    | aggregate_results

  emit:
  output_ch
}

workflow add_methods {
  take: input_ch
  main:
  output_ch = Channel.fromList(methods.keySet())
    | combine(input_ch)

    // generate combined id for method_id and dataset_id
    | pmap{method_id, dataset_id, data ->
      def new_id = dataset_id + "." + method_id
      def new_data = data.clone() + [method_id: method_id]
      new_data.remove("id")
      [new_id, new_data]
    }
  emit: output_ch
}

workflow check_filtered_normalization_id {
  take: input_ch
  main:
  output_ch = input_ch
    | pfilter{id, data ->
      data = data.clone()
      def method = methods[data.method_id]
      def preferred = method.config.functionality.info.preferred_normalization
      // if a method is just using the counts, we can use any normalization method
      if (preferred == "counts") {
        preferred = "log_cpm"
      }
      data.normalization_id == preferred
    }
  emit: output_ch
}

workflow controls_can_cheat {
  take: input_ch
  main:
  output_ch = input_ch
    | pmap{id, data, passthrough ->
      def method = methods[data.method_id]
      def method_type = method.config.functionality.info.method_type
      def new_data = data.clone()
      // if (method_type != "method") {
      //   new_data = new_data + [input_solution: passthrough.metric.input_solution]
      // }
      [id, new_data, passthrough]
    }
  emit: output_ch
}

workflow run_methods {
  take: input_ch
  main:
    // generate one channel per method
    method_chs = methods.collect { method_id, method_module ->
        input_ch
          | filter{it[1].method_id == method_id}
          | method_module
      }
    // mix all results
    output_ch = method_chs[0].mix(*method_chs.drop(1))

  emit: output_ch
}

workflow run_metrics {
  take: input_ch
  main:

  output_ch = input_ch
    | (rmse & trustworthiness & density)
    | mix

  emit: output_ch
}

workflow aggregate_results {
  take: input_ch
  main:

  output_ch = input_ch
    | toSortedList
    | filter{ it.size() > 0 }
    | map{ it -> 
      [ "combined", it.collect{ it[1] } ] + it[0].drop(2) 
    }
    | getWorkflowArguments(key: "output")
    | extract_scores.run(
        auto: [ publish: true ]
    )

  emit: output_ch
}