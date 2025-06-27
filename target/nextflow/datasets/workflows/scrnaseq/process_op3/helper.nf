Map findArgumentSchema(Map config, String argument_id) {
  def argument_groups =
    (config.argument_groups ?: []) +
    [
      arguments: config.arguments ?: []
    ]

  def schema_value = argument_groups.findResult{ gr ->
    gr.arguments.find { arg ->
      arg.name == ("--" + argument_id)
    }
  }
  return schema_value
}

Boolean checkItemAllowed(String item, List include, List exclude, String includeArgName, String excludeArgName) {

  // Throw an error if both include and exclude lists are provided
  if (include != null && exclude != null) {
      throw new Exception("Cannot define both ${includeArgName} and ${excludeArgName}")
  }

  if (include) {
    return include.contains(item)
  }
  if (exclude) {
    return !exclude.contains(item)
  }

  return true
}
