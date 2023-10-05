// add custom tracer to nextflow to capture exit codes, memory usage, cpu usage, etc.
traces = initializeTracer()

workflow run_wf {
  take:
  input_ch

  main:

  normalization_methods = [log_cp, sqrt_cp, l1_sqrt, log_scran_pooling]


  dataset_ch = input_ch

    // fetch data from legacy openproblems
    | openproblems_v1_multimodal.run(
      fromState: { id, state -> 
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
          dataset_organism: state.dataset_organism
        ]
      },
      toState: [ 
        raw_mod1: "output_mod1",
        raw_mod2: "output_mod2"
      ]
    )

  sampled_dataset_ch = dataset_ch
    | filter{ id, state -> state.do_subsample }
    | subsample.run(
      fromState: { id, state ->
        [
          input: state.raw_mod1,
          input_mod2: state.raw_mod2,
          n_obs: state.n_obs,
          n_vars: state.n_vars,
          keep_features: state.keep_features,
          keep_celltype_categories: state.keep_celltype_categories,
          keep_batch_categories: state.keep_batch_categories,
          even: state.even,
          seed: state.seed,
          output_mod2: '$id.$key.output_mod2.h5ad' // set value for optional output
        ]
      },
      toState: [
        raw_mod1: "output",
        raw_mod2: "output_mod2"
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
          input: state.raw_mod1,
          output: '$id.$key.output_mod1.h5ad'
        ]
      },
      toState: { id, output, state, config -> 
        state + [
          normalization_id: config.functionality.name,
          normalized_mod1: output.output
        ]
      }
    )
    // run normalization methods on second modality
    | runComponents(
      components: normalization_methods,
      filter: { id, state, config ->
        config.functionality.name == state.normalization_id
      },
      fromState: { id, state, config ->
        [
          input: state.raw_mod2,
          output: '$id.$key.output_mod2.h5ad'
        ]
      },
      toState: [normalized_mod2: "output"]
    )

    | svd.run(
      fromState: { id, state ->
        [
          input: state.normalized_mod1,
          input_mod2: state.normalized_mod2,
          output: '$id.$key.output_mod1.h5ad',
          output_mod2: '$id.$key.output_mod2.h5ad'
        ]
      },
      toState: [
        svd_mod1: "output",
        svd_mod2: "output_mod2"
      ]
    )

    | hvg.run(
      fromState: [ input: "svd_mod1" ],
      toState: [ hvg_mod1: "output" ]
    )

    | hvg.run(
      fromState: [ input: "svd_mod2" ],
      toState: [ hvg_mod2: "output" ]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.hvg_mod1,
          meta: state.output_meta_mod1 ?: '$id.$key.output_meta_mod1.yaml',
          output: state.output_dataset_mod1 ?: '$id.$key.output_dataset_mod1.h5ad',
          checks: null
        ]
      },
      toState: [ dataset_mod1: "output", meta_mod1: "meta" ]
    )

    | check_dataset_schema.run(
      fromState: { id, state ->
        [
          input: state.hvg_mod2,
          meta: state.output_meta_mod2 ?: '$id.$key.output_meta_mod2.yaml',
          output: state.output_dataset_mod2 ?: '$id.$key.output_dataset_mod2.h5ad',
          checks: null
        ]
      },
      toState: [ dataset_mod2: "output", meta_mod2: "meta" ]
    )

    // only output the files for which an output file was specified
    | setState{ id, state ->
      [
        "output_dataset_mod1" : state.output_dataset_mod1 ? state.dataset_mod1: null,
        "output_dataset_mod2" : state.output_dataset_mod2 ? state.dataset_mod2: null,
        "output_meta_mod1" : state.output_meta_mod1 ? state.meta_mod1: null,
        "output_meta_mod2" : state.output_meta_mod2 ? state.meta_mod2: null
      ]
    }

  emit:
  output_ch
}

// store the trace log in the publish dir
workflow.onComplete {
  def publish_dir = getPublishDir()

  writeJson(traces, file("$publish_dir/traces.json"))
}