nextflow.enable.dsl=2

sourceDir = params.rootDir + "/src"
targetDir = params.rootDir + "/target/nextflow"

// dataset loaders
include { openproblems_v1 } from "$targetDir/datasets/loaders/openproblems_v1/main.nf"

// normalization methods
include { log_cp } from "$targetDir/datasets/normalization/log_cp/main.nf"
include { log_scran_pooling } from "$targetDir/datasets/normalization/log_scran_pooling/main.nf"
include { sqrt_cp } from "$targetDir/datasets/normalization/sqrt_cp/main.nf"
include { l1_sqrt } from "$targetDir/datasets/normalization/l1_sqrt/main.nf"

// dataset processors
include { subsample } from "$targetDir/datasets/processors/subsample/main.nf"
include { pca } from "$targetDir/datasets/processors/pca/main.nf"
include { hvg } from "$targetDir/datasets/processors/hvg/main.nf"
include { knn } from "$targetDir/datasets/processors/knn/main.nf"
include { check_dataset_schema } from "$targetDir/common/check_dataset_schema/main.nf"

// helper functions
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from sourceDir + "/wf_utils/WorkflowHelper.nf"
include { publishState; runComponents; joinStates; initializeTracer; writeJson; getPublishDir; setState } from sourceDir + "/wf_utils/BenchmarkHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

normalization_methods = [log_cp, sqrt_cp, l1_sqrt, log_scran_pooling]

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
    | publishState([:])
}

workflow run_wf {
  take:
  input_ch

  main:
  dataset_ch = input_ch
    | preprocessInputs(config: config)

    // fetch data from legacy openproblems
    | openproblems_v1.run(
      fromState: { id, state -> 
        def output_filename =
          (!state.do_subsample && state.output_raw) ? 
          state.output_raw :
          '$id.$key.output_raw.h5ad'
        [
          dataset_id: state.dataset_id,
          obs_celltype: state.obs_celltype,
          obs_batch: state.obs_batch,
          obs_tissue: state.obs_tussue,
          layer_counts: state.layer_counts,
          sparse: state.sparse,
          dataset_name: state.dataset_name,
          data_url: state.data_url,
          data_reference: state.data_reference, 
          dataset_summary: state.dataset_summary,
          dataset_description: state.dataset_description,
          dataset_organism: state.dataset_organism,
          output: output_filename
        ]
      },
      toState: [ 
        raw: "output"
      ]
    )

  sampled_dataset_ch = dataset_ch
    | filter{ id, state -> state.do_subsample }
    | subsample.run(
      fromState: { id, state ->
        [
          input: state.raw,
          n_obs: state.n_obs,
          n_vars: state.n_vars,
          keep_features: state.keep_features,
          keep_celltype_categories: state.keep_celltype_categories,
          keep_batch_categories: state.keep_batch_categories,
          even: state.even,
          seed: state.seed,
          output: state.output_raw ?: '$id.$key.output_raw.h5ad',
          output_mod2: null
        ]
      },
      toState: [
        raw: "output"
      ]
    )
  notsampled_dataset_ch = dataset_ch
    | filter{ id, state -> !state.do_subsample }
  
  output_ch = sampled_dataset_ch
    | concat(notsampled_dataset_ch)

    // run normalization methods
    | runComponents(
      components: normalization_methods,
      id: { id, state, config ->
        if (state.normalization_methods.size() > 1) {
          id + "/" + config.functionality.name
        } else {
          id
        }
      },
      filter: { id, state, config ->
        config.functionality.name in state.normalization_methods
      },
      fromState: { id, state, config -> 
        [
          input: state.raw,
          output: state.output_normalized ?: '$id.$key.output_normalized.h5ad'
        ]
      },
      toState: { id, output, state, config -> 
        state + [
          normalization_id: config.functionality.name,
          normalized: output.output
        ]
      }
    )

    | pca.run(
      fromState: { id, state ->
        [
          input: state.normalized,
          output: state.output_pca ?: '$id.$key.output_pca.h5ad'
        ]
      },
      toState: [ pca: "output" ]
    )

    | hvg.run(
      fromState: { id, state -> 
        [
          input: state.pca,
          output: state.output_hvg ?: '$id.$key.output_hvg.h5ad'
        ]
      },
      toState: [ hvg: "output" ]
    )

    | knn.run(
      fromState: { id, state ->
        [
          input: state.hvg,
          output: state.output_knn ?: '$id.$key.output_knn.h5ad'
        ]
      },
      toState: [ knn: "output" ]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.knn,
          meta: state.output_meta ?: '$id.$key.output_meta.yaml',
          output: state.output_dataset ?: '$id.$key.output_dataset.h5ad',
          checks: null
        ]
      },
      toState: [ dataset: "output", meta: "meta" ]
    )

    // only output the files for which an output file was specified
    | setState{ id, state ->
      [
        "output_dataset": state.output_dataset ? state.dataset : null,
        "output_meta": state.output_meta ? state.meta : null,
        "output_raw": state.output_raw ? state.raw : null,
        "output_normalized": state.output_normalized ? state.normalized : null,
        "output_pca": state.output_pca ? state.pca : null,
        "output_hvg": state.output_hvg ? state.hvg : null,
        "output_knn": state.output_knn ? state.knn : null
      ]
    }

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
  // writeJson(normalization_methods.collect{it.config}, file("$publish_dir/normalization_methods.json"))
}