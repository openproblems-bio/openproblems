workflow run_wf {
  take: input_ch
  main:
  output_ch = input_ch
    | novel_train.run(
      fromState: ["input_train_mod1", "input_train_mod2"],
      toState: ["input_model": "output", "input_transform": "output_transform", "output_train_mod2": "output_train_mod2"]
    )
    | novel_predict.run(
      fromState: { id, state ->
        [
          "input_train_mod2": state.output_train_mod2,
          "input_test_mod1": state.input_test_mod1,
          "input_model": state.input_model, 
          "input_transform": state.input_transform,
          "output": state.output]},
      toState: ["output": "output"]
    )

    | map { tup ->
      [tup[0], [output: tup[1].output]]
    }

  emit: output_ch
}