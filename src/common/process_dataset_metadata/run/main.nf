workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | yaml_to_json.run(
      fromState: ["input"],
      toState: ["output"]
    )

    | setState(["output"])

    emit:
    output_ch
}