nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import control methods
include { no_denoising } from "$targetDir/denoising/control_methods/no_denoising/main.nf"
include { perfect_denoising } from "$targetDir/denoising/control_methods/perfect_denoising/main.nf"

// import methods
include { alra } from "$targetDir/denoising/methods/alra/main.nf"
include { dca } from "$targetDir/denoising/methods/dca/main.nf"
include { knn_smoothing } from "$targetDir/denoising/methods/knn_smoothing/main.nf"
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

// construct a map of methods (id -> method_module)
methods = [ no_denoising, perfect_denoising, alra, dca, knn_smoothing, magic]
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
    | setWorkflowArguments(
      method: ["input_train", "input_test"],
      metric: ["input_test"],
      output: ["output"]
    )

    // multiply events by the number of method
    | add_methods

    // add input_solution to data for the positive controls
    | controls_can_cheat

    // run methods
    | getWorkflowArguments(key: "method")
    | run_methods

    // construct tuples for metrics
    | pmap{ id, file, passthrough ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      // combine prediction with solution
      def newData = [ input_denoised: file, input_test: passthrough.metric.input_test ]
      [ newId, newData, passthrough ]
    }
    
    // run metrics
    | getWorkflowArguments(key: "metric")
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

// workflow check_filtered_normalization_id {
//   take: input_ch
//   main:
//   output_ch = input_ch
//     | pfilter{id, data ->
//       data = data.clone()
//       def method = methods[data.method_id]
//       def preferred = method.config.functionality.info.preferred_normalization
//       // if a method is just using the counts, we can use any normalization method
//       if (preferred == "counts") {
//         preferred = "log_cpm"
//       }
//       data.normalization_id == preferred
//     }
//   emit: output_ch
// }

workflow controls_can_cheat {
  take: input_ch
  main:
  output_ch = input_ch
    | pmap{id, data, passthrough ->
      def method = methods[data.method_id]
      def method_type = method.config.functionality.info.method_type
      def new_data = data.clone()
      if (method_type != "method") {
        new_data = new_data + [input_test: passthrough.metric.input_test]
      }
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
    | (mse & poisson)
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