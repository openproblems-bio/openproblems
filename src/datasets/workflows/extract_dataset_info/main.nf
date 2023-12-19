workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // extract the dataset metadata
    | check_dataset_schema.run(
      fromState: [input: "input"],
      toState: { id, output, state ->
        def dataset_uns = (new org.yaml.snakeyaml.Yaml().load(output.meta)).uns
        state + [dataset_uns: dataset_uns]
      }
    )

    // only keep one of the normalization methods
    | filter{ id, state ->
      if (state.filter_normalization_id) {
        state.filter_normalization_id.contains(state.dataset_uns.normalization_id)
      } else {
        true
      }
    }

    | joinStates { ids, states ->
      // remove normalization id
      // TODO: make this optional through a parameter?
      def dataset_uns = states.collect{state ->
        def uns = state.dataset_uns.clone()
        uns.remove("normalization_id")
        uns
      }

      // store data as yaml
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      def new_state = [
        output: dataset_uns_file, 
        _meta: [join_id: ids[0]]
      ]
      ["output", new_state]
    }


  emit:
  output_ch
}
