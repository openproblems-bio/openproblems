Map findArgumentSchema(Map config, String argument_id) {
  def argument_groups =
    (config.functionality.argument_groups ?: []) +
    [
      arguments: config.functionality.arguments ?: []
    ]

  def schema_value = argument_groups.findResult{ gr ->
    gr.arguments.find { arg ->
      arg.name == ("--" + argument_id)
    }
  }
  return schema_value
}
