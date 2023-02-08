nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import methods
include { bbknn } from "$targetDir/batch_integration/graph/methods/bbknn/main.nf"
include { combat } from "$targetDir/batch_integration/graph/methods/combat/main.nf"
include { scanorama_embed } from "$targetDir/batch_integration/graph/methods/scanorama_embed/main.nf"
include { scanorama_feature } from "$targetDir/batch_integration/graph/methods/scanorama_feature/main.nf"
include { scvi } from "$targetDir/batch_integration/graph/methods/scvi/main.nf"

// import metrics
include { clustering_overlap } from "$targetDir/batch_integration/graph/metrics/clustering_overlap/main.nf"

// tsv generation component
include { extract_scores } from "$targetDir/common/extract_scores/main.nf"

// import helper functions
include { readConfig; viashChannel; helpMessage } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from sourceDir + "/wf_utils/DataflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

      // split params for downstream components
    | setWorkflowArguments(
      method: ["input"],
      metric: [],
      output: ["output"]
    )

    // multiply events by the number of method
    | add_methods

    // run methods
    | getWorkflowArguments(key: "method")
    | run_methods

    // construct tuples for metrics
    | pmap{ id, file, passthrough ->
      // derive unique ids from output filenames
      def newId = file.getName().replaceAll(".output.*", "")
      // combine prediction with solution
      def newData = [ input: file ]
      [ newId, newData, passthrough ]
    }
    
    // run metrics
    | getWorkflowArguments(key: "metric")
    | run_metrics
    
    // convert to tsv  
    | aggregate_results

  emit:
  output_ch
/*
        | (bbknn & combat & scvi & scanorama_embed & scanorama_feature)
        | mix
        | toSortedList
        | view { "toSortedList $it" }
/*
        | map{ it -> [ "combined", [ input: it.collect{ it[1] } ] ] }
        | (ari & nmi)
        | extract_scores.run(
            auto: [ publish: true ]
        )
*/
}

/*******************************************************
*                    Sub workflows                     *
*******************************************************/

// construct a map of methods (id -> method_module)
methods = [ bbknn, combat, scanorama_embed, scanorama_feature, scvi]
  .collectEntries{method ->
    [method.config.functionality.name, method]
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
    | (ari & nmi)
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
