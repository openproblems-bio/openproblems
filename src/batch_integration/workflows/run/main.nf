nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// import methods
include { bbknn } from "$targetDir/batch_integration/methods_graph/bbknn/main.nf"
include { combat } from "$targetDir/batch_integration/methods_feature/combat/main.nf"
include { scanorama_embed } from "$targetDir/batch_integration/methods_embedding/scanorama_embed/main.nf"
include { scanorama_feature } from "$targetDir/batch_integration/methods_feature/scanorama_feature/main.nf"
include { scvi } from "$targetDir/batch_integration/methods_graph/scvi/main.nf"

// import transformers
include ( feature_to_embed ) from "$targetDir/batch_integration/methods_feature/feature_to_embed/main.nf"
include ( embed_to_graph ) from "$targetDir/batch_integration/methods_graph/embed_to_graph/main.nf"

// import metrics
include { clustering_overlap } from "$targetDir/batch_integration/metrics_graph/clustering_overlap/main.nf"
include { asw_batch } from "$targetDir/batch_integration/metrics_embedding/asw_batch/main.nf"
include { asw_label } from "$targetDir/batch_integration/metrics_embedding/asw_label/main.nf"
include { ccc } from "$targetDir/batch_integration/metrics_embedding/cel_cycle_conservation/main.nf"
include { pcr } from "$targetDir/batch_integration/metrics_embedding/pcr/main.nf"

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

  // run feature methods
  meth_feature = input_ch
    | (combat & scanorama_feature)
    | mix
  // run embed methods
  meth_embed = input_ch
    | (scanorama_embed & scvi)
    | mix
  // run graph methods
  meth_graph = input_ch
    | (bbknn)

  // apply feature metrics on feature outputs
  // metr_feat = meth_feature
  //   | (asw_batch & asw_label & cell_cycle_conservation & pcr)

  // convert feature outputs to embedding outputs
  meth_feat_to_embed = meth_feature
    | feature_to_embed

  // apply embedding metrics to embedding outputs
  metr_embed = meth_embed
    | mix(meth_feat_to_embed)
    | (asw_batch & asw_label & cell_cycle_conservation & pcr)
    | mix
  
  // convert embedding outputs to graph outputs
  meth_embed_to_graph = meth_embed
    | mix(meth_feat_to_embed)
    | embed_to_graph
  
  // apply graph metrics to graph outputs
  metr_graph = meth_graph
    | mix(meth_embed_to_graph)
    | (clustering_overlap)

  
  output_ch = metr_embed.mix(metr_graph)
    
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
