////////////////////////////
// VDSL3 helper functions //
////////////////////////////

// helper file: 'src/main/resources/io/viash/platforms/nextflow/arguments/_processArgumentGroup.nf'
def _processArgumentGroup(argumentGroups, name, arguments) {
  def argNamesInGroups = argumentGroups.collectMany{it.arguments.findAll{it instanceof String}}.toSet()

  // Check if 'arguments' is in 'argumentGroups'. 
  def argumentsNotInGroup = arguments.findAll{arg -> !(argNamesInGroups.contains(arg.plainName))}

  // Check whether an argument group of 'name' exists.
  def existing = argumentGroups.find{gr -> name == gr.name}

  // if there are no arguments missing from the argument group, just return the existing group (if any)
  if (argumentsNotInGroup.isEmpty()) {
    return existing == null ? [] : [existing]
  
  // if there are missing arguments and there is an existing group, add the missing arguments to it
  } else if (existing != null) {
    def newEx = existing.clone()
    newEx.arguments.addAll(argumentsNotInGroup.findAll{it !instanceof String})
    return [newEx]

  // else create a new group
  } else {
    def newEx = [name: name, arguments: argumentsNotInGroup.findAll{it !instanceof String}]
    return [newEx]
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/arguments/_processArgument.nf'
def _processArgument(arg) {
  arg.multiple = arg.multiple != null ? arg.multiple : false
  arg.required = arg.required != null ? arg.required : false
  arg.direction = arg.direction != null ? arg.direction : "input"
  arg.multiple_sep = arg.multiple_sep != null ? arg.multiple_sep : ":"
  arg.plainName = arg.name.replaceAll("^-*", "")

  if (arg.type == "file") {
    arg.must_exist = arg.must_exist != null ? arg.must_exist : true
    arg.create_parent = arg.create_parent != null ? arg.create_parent : true
  }

  // add default values to output files which haven't already got a default
  if (arg.type == "file" && arg.direction == "output" && arg.default == null) {
    def mult = arg.multiple ? "_*" : ""
    def extSearch = ""
    if (arg.default != null) {
      extSearch = arg.default
    } else if (arg.example != null) {
      extSearch = arg.example
    }
    if (extSearch instanceof List) {
      extSearch = extSearch[0]
    }
    def extSearchResult = extSearch.find("\\.[^\\.]+\$")
    def ext = extSearchResult != null ? extSearchResult : ""
    arg.default = "\$id.\$key.${arg.plainName}${mult}${ext}"
    if (arg.multiple) {
      arg.default = [arg.default]
    }
  }

  if (!arg.multiple) {
    if (arg.default != null && arg.default instanceof List) {
      arg.default = arg.default[0]
    }
    if (arg.example != null && arg.example instanceof List) {
      arg.example = arg.example[0]
    }
  }

  if (arg.type == "boolean_true") {
    arg.default = false
  }
  if (arg.type == "boolean_false") {
    arg.default = true
  }

  arg
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/arguments/processInputsOutputs.nf'
def typeCheck(String stage, Map par, Object value, String id, String key) {
  // expectedClass will only be != null if value is not of the expected type
  def expectedClass = null
  
  if (!par.required && value == null) {
    expectedClass = null
  } else if (par.multiple) {
    if (value instanceof List) {
      try {
        value = value.collect { listVal ->
          typeCheck(stage, par + [multiple: false], listVal, id, key)
        }
      } catch (Exception e) {
        expectedClass = "List[${par.type}]"
      }
    } else {
      expectedClass = "List[${par.type}]"
    }
  } else if (par.type == "string") {
    if (value instanceof GString) {
      value = value.toString()
    }
    expectedClass = value instanceof String ? null : "String"
  } else if (par.type == "integer") {
    expectedClass = value instanceof Integer ? null : "Integer"
  } else if (par.type == "long") {
    if (value instanceof Integer) {
      value = value.toLong()
    }
    expectedClass = value instanceof Long ? null : "Long"
  } else if (par.type == "boolean") {
    expectedClass = value instanceof Boolean ? null : "Boolean"
  } else if (par.type == "file") {
    if (stage == "output" || par.direction == "input") {
      if (value instanceof File) {
        value = value.toPath()
      }
      expectedClass = value instanceof Path ? null : "Path"
    } else { // stage == "input" && par.direction == "output"
      if (value instanceof GString) {
        value = value.toString()
      }
      expectedClass = value instanceof String ? null : "String"
    }
  } else {
    expectedClass = par.type
  }

  if (expectedClass != null) {
    error "Error in module '${key}' id '${id}': ${stage} argument '${par.plainName}' has the wrong type. " +
      "Expected type: ${expectedClass}. Found type: ${value.getClass()}"
  }
  
  return value
}

Map processInputs(Map inputs, Map config, String id, String key) {
  if (!workflow.stubRun) {
    config.functionality.allArguments.each { arg ->
      if (arg.required) {
        assert inputs.containsKey(arg.plainName) && inputs.get(arg.plainName) != null : 
          "Error in module '${key}' id '${id}': required input argument '${arg.plainName}' is missing"
      }
    }

    inputs = inputs.collectEntries { name, value ->
      def par = config.functionality.allArguments.find { it.plainName == name && (it.direction == "input" || it.type == "file") }
      assert par != null : "Error in module '${key}' id '${id}': '${name}' is not a valid input argument"

      value = typeCheck("input", par, value, id, key)

      [ name, value ]
    }
  }
  return inputs
}

Map processOutputs(Map outputs, Map config, String id, String key) {
  if (!workflow.stubRun) {
    config.functionality.allArguments.each { arg ->
      if (arg.direction == "output" && arg.required) {
        assert outputs.containsKey(arg.plainName) && outputs.get(arg.plainName) != null : 
          "Error in module '${key}' id '${id}': required output argument '${arg.plainName}' is missing"
      }
    }

    outputs = outputs.collectEntries { name, value ->
      def par = config.functionality.allArguments.find { it.plainName == name && it.direction == "output" }
      assert par != null : "Error in module '${key}' id '${id}': '${name}' is not a valid output argument"
      
      value = typeCheck("output", par, value, id, key)
      
      [ name, value ]
    }
  }
  return outputs
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/applyConfig.nf'

/**
 * Apply the argument settings specified in a Viash config to a list of parameter sets.
 *    - Split the parameter values according to their seperator if 
 *       the parameter accepts multiple values
 *    - Cast the parameters to their corect types.
 *    - Assertions:
 *        ~ Check if any unknown parameters are found
 *        ~ Check if the ID of the parameter set is unique across all sets.
 * 
 * @return The input parameters that have been processed.
 */

List<Tuple> applyConfig(List<Tuple> parameterSets, Map config){
  def processedparameterSets = parameterSets.collect({ parameterSet ->
    def id = parameterSet[0]
    def paramValues = parameterSet[1]
    def passthrough = parameterSet.drop(2)
    def processedSet = applyConfigToOneParameterSet(paramValues, config)
    [id, processedSet] + passthrough
  })

  _checkUniqueIds(processedparameterSets)
  return processedparameterSets
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/applyConfigToOneParameterSet.nf'
/**
 * Cast parameters to the correct type as defined in the Viash config
 *
 * @param parValues A Map of input arguments.
 *
 * @return The input arguments that have been cast to the type from the viash config.
 */
private Map<String, Object> _castParamTypes(Map<String, Object> parValues, Map config) {
  // Cast the input to the correct type according to viash config
  def castParValues = parValues.collectEntries({ parName, parValue ->
    def paramSettings = config.functionality.allArguments.find({it.plainName == parName})
    // dont parse parameters like publish_dir ( in which case paramSettings = null)
    def parType = paramSettings ? paramSettings.get("type", null) : null

    // turn parValue into a list, if it isn't one already
    if (parValue !instanceof Collection) {
      parValue = [parValue]
    }

    if (parType == "file" && paramSettings.get("direction", "input") == "input") {
      parValue = parValue.collectMany{ parValueItem ->
        if (parValueItem instanceof String) {
          parValueItem = file(parValueItem)
        }
        if (parValueItem !instanceof Collection) {
          parValueItem = [parValueItem]
        }
        parValueItem
      }
      // cast Paths to File? Or vice versa?
      // parValue = parValue.collect{
      //   if (path instanceof Path) {
      //     [path.toFile()]
      //   } else {
      //     [path]
      //   }
      // }
    } else if (parType == "integer") {
      parValue = parValue.collect{it as Integer}
    } else if (parType == "double") {
      parValue = parValue.collect{it as Double}
    } else if (parType == "boolean" || 
                parType == "boolean_true" || 
                parType == "boolean_false") {
      parValue = parValue.collect{it as Boolean}
    }

    // simplify list to value if need be
    if (paramSettings && !paramSettings.multiple) {
      assert parValue.size() == 1 : 
        "Error: argument ${parName} has too many values.\n" +
        "  Expected amount: 1. Found: ${parValue.size()}"
      parValue = parValue[0]
    }

    [parName, parValue]
  })
  return castParValues
}

/**
 * Apply the argument settings specified in a Viash config to a single parameter set.
 *    - Split the parameter values according to their seperator if 
 *       the parameter accepts multiple values
 *    - Cast the parameters to their corect types.
 *    - Assertions:
 *        ~ Check if any unknown parameters are found
 * 
 * @param paramValues A Map of parameter to be processed. All parameters must 
 *                    also be specified in the Viash config.
 * @param config: A Map of the Viash configuration. This Map can be generated from 
 *                the config file using the readConfig() function.
 * @return The input parameters that have been processed.
 */
Map<String, Object> applyConfigToOneParameterSet(Map<String, Object> paramValues, Map config){
  def splitParamValues = _splitParams(paramValues, config)
  def castParamValues = _castParamTypes(splitParamValues, config)

  // Check if any unexpected arguments were passed
  def knownParams = config.functionality.allArguments.collect({it.plainName}) + ["publishDir", "publish_dir"]
  castParamValues.each({parName, parValue ->
      assert parName in knownParams: "Unknown parameter. Parameter $parName should be in $knownParams"
  })
  return castParamValues
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/channelFromParams.nf'
/**
 * Resolve the file paths in the parameters relative to given path
 *
 * @param paramList A Map containing parameters to process.
 *                  This function assumes that files are still of type String.
 * @param config A Map of the Viash configuration. This Map can be generated from the config file
 *               using the readConfig() function.
 * @param relativeTo path of a file to resolve the parameters values to.
 *
 * @return A map of parameters where the location of the input file parameters have been resolved
 *         resolved relatively to the provided path.
 */
private Map<String, Object> _resolvePathsRelativeTo(Map paramList, Map config, String relativeTo) {
  paramList.collectEntries { parName, parValue ->
    argSettings = config.functionality.allArguments.find{it.plainName == parName}
    if (argSettings && argSettings.type == "file" && argSettings.direction == "input") {
      if (parValue instanceof Collection) {
        parValue = parValue.collect({path -> 
          path !instanceof String ? path : file(_getChild(relativeTo, path))
        })
      } else {
        parValue = parValue !instanceof String ? path : file(_getChild(relativeTo, parValue))
      }
    }
    [parName, parValue]
  }
}

/**
 * Parse nextflow parameters based on settings defined in a viash config 
 * and return a nextflow channel.
 *
 * @param params Input parameters from nextflow.
 * @param config A Map of the Viash configuration. This Map can be generated from the config file
 *               using the readConfig() function.
 *
 * @return A list of parameter sets that were parsed from the 'param_list' argument value.
 */
private List<Tuple2<String, Map>> _parseParamListArguments(Map params, Map config){
  // first try to guess the format (if not set in params)
  def paramListFormat = _guessParamListFormat(params)

  // get the correct parser function for the detected params_list format
  def paramListParsers = [ 
    "csv": {[it, readCsv(it)]},
    "json": {[it, readJson(it)]},
    "yaml": {[it, readYaml(it)]},
    "yaml_blob": {[null, readYamlBlob(it)]},
    "asis": {[null, it]},
    "none": {[null, [[:]]]}
  ]
  assert paramListParsers.containsKey(paramListFormat):
    "Format of provided --param_list not recognised.\n" +
    "You can use '--param_list_format' to manually specify the format.\n" +
    "Found: '$paramListFormat'. Expected: one of 'csv', 'json', "+
    "'yaml', 'yaml_blob', 'asis' or 'none'"
  def paramListParser = paramListParsers.get(paramListFormat)

  // fetch multi param inputs
  def paramListOut = paramListParser(params.containsKey("param_list") ? params.param_list : "")
  // multiFile is null if the value passed to param_list was not a file (e.g a blob)
  // If the value was indeed a file, multiFile contains the location that file (used later).
  def paramListFile = paramListOut[0]
  def paramSets = paramListOut[1] // these are the actual parameters from reading the blob/file

  // data checks
  assert paramSets instanceof List: "--param_list should contain a list of maps"
  for (value in paramSets) {
    assert value instanceof Map: "--param_list should contain a list of maps"
  }

  // id is argument
  def idIsArgument = config.functionality.allArguments.find({it.plainName == "id"}) != null

  // Reformat from List<Map> to List<Tuple2<String, Map>> by adding the ID as first element of a Tuple2
  paramSets = paramSets.collect({ paramValues ->
    def paramId = paramValues.id
    if (!idIsArgument) {
      paramValues = paramValues.findAll{k, v -> k != "id"}
    }
    [paramId, paramValues]
  })

  // Split parameters with 'multiple: true'
  paramSets = paramSets.collect({ id, paramValues ->
    def splitParamValues = _splitParams(paramValues, config)
    [id, splitParamValues]
  })
  
  // The paths of input files inside a param_list file may have been specified relatively to the
  // location of the param_list file. These paths must be made absolute.
  if (paramListFile){
    paramSets = paramSets.collect({ id, paramValues ->
      def relativeParamValues = _resolvePathsRelativeTo(paramValues, config, paramListFile)
      [id, relativeParamValues]
    })
  }

  return paramSets
}

/**
 * Parse nextflow parameters based on settings defined in a viash config.
 * Return a list of parameter sets, each parameter set corresponding to 
 * an event in a nextflow channel. The output from this function can be used
 * with Channel.fromList to create a nextflow channel with Vdsl3 formatted 
 * events.
 *
 * This function performs:
 *   - A filtering of the params which can be found in the config file.
 *   - Process the params_list argument which allows a user to to initialise 
 *     a Vsdl3 channel with multiple parameter sets. Possible formats are 
 *     csv, json, yaml, or simply a yaml_blob. A csv should have column names 
 *     which correspond to the different arguments of this pipeline. A json or a yaml
 *     file should be a list of maps, each of which has keys corresponding to the
 *     arguments of the pipeline. A yaml blob can also be passed directly as a parameter.
 *     When passing a csv, json or yaml, relative path names are relativized to the
 *     location of the parameter file.
 *   - Combine the parameter sets into a vdsl3 Channel.
 *
 * @param params Input parameters. Can optionaly contain a 'param_list' key that
 *               provides a list of arguments that can be split up into multiple events
 *               in the output channel possible formats of param_lists are: a csv file, 
 *               json file, a yaml file or a yaml blob. Each parameters set (event) must
 *               have a unique ID.
 * @param config A Map of the Viash configuration. This Map can be generated from the config file
 *               using the readConfig() function.
 * 
 * @return A list of parameters with the first element of the event being
 *         the event ID and the second element containing a map of the parsed parameters.
 */
 
private List<Tuple2<String, Map<String, Object>>> _paramsToParamSets(Map params, Map config){
  /* parse regular parameters (not in param_list)  */
  /*************************************************/
  def globalParams = config.functionality.allArguments
    .findAll { params.containsKey(it.plainName) }
    .collectEntries { [ it.plainName, params[it.plainName] ] }
  def globalID = params.get("id", null)
  def globalParamsValues = applyConfigToOneParameterSet(globalParams, config)

  /* process params_list arguments */
  /*********************************/
  def paramSets = _parseParamListArguments(params, config)
  def parameterSetsWithConfigApplied = applyConfig(paramSets, config)

  /* combine arguments into channel */
  /**********************************/
  def processedParams = parameterSetsWithConfigApplied.indexed().collect{ index, paramSet ->
    def id = paramSet[0]
    def parValues = paramSet[1]
    id = [id, globalID].find({it != null}) // first non-null element
  
    if (workflow.stubRun) {
      // if stub run, explicitly add an id if missing
      id = id ? id : "stub" + index
    }
    assert id != null: "Each parameter set should have at least an 'id'"
    // Add regular parameters together with parameters passed with 'param_list'
    def combinedArgsValues = globalParamsValues + parValues

    // Remove parameters which are null, if the default is also null
    combinedArgsValues = combinedArgsValues.collectEntries{paramName, paramValue ->
      parameterSettings = config.functionality.allArguments.find({it.plainName == paramName})
      if ( paramValue != null || parameterSettings.get("default", null) != null ) {
        [paramName, paramValue]
      }
    }
    [id, combinedArgsValues]
  }

  // Check if ids (first element of each list) is unique
  _checkUniqueIds(processedParams)
  return processedParams
}

/**
 * Parse nextflow parameters based on settings defined in a viash config 
 * and return a nextflow channel.
 * 
 * @param params Input parameters. Can optionaly contain a 'param_list' key that
 *               provides a list of arguments that can be split up into multiple events
 *               in the output channel possible formats of param_lists are: a csv file, 
 *               json file, a yaml file or a yaml blob. Each parameters set (event) must
 *               have a unique ID.
 * @param config A Map of the Viash configuration. This Map can be generated from the config file
 *               using the readConfig() function.
 * 
 * @return A nextflow Channel with events. Events are formatted as a tuple that contains 
 *         first contains the ID of the event and as second element holds a parameter map.
 *       
 *
 */
def channelFromParams(Map params, Map config) {
  processedParams = _paramsToParamSets(params, config)
  return Channel.fromList(processedParams)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/_checkUniqueIds.nf'

/**
 * Check if the ids are unique across parameter sets
 *
 * @param parameterSets a list of parameter sets.
 */
private void _checkUniqueIds(List<Tuple2<String, Map<String, Object>>> parameterSets) {
  def ppIds = parameterSets.collect{it[0]}
  assert ppIds.size() == ppIds.unique().size() : "All argument sets should have unique ids. Detected ids: $ppIds"
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/checkUniqueIds.nf'
def checkUniqueIds(Map args) {
  def stopOnError = args.stopOnError == null ? args.stopOnError : true

  def idChecker = new IDChecker()

  return filter { tup ->
    if (!idChecker.observe(tup[0])) {
      if (stopOnError) {
        error "Duplicate id: ${tup[0]}"
      } else {
        log.warn "Duplicate id: ${tup[0]}, removing duplicate entry"
        return false
      }
    }
    return true
  }
}
// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/_getChild.nf'

// helper functions for reading params from file //
def _getChild(parent, child) {
  if (child.contains("://") || java.nio.file.Paths.get(child).isAbsolute()) {
    child
  } else {
    def parentAbsolute = java.nio.file.Paths.get(parent).toAbsolutePath().toString()
    parentAbsolute.replaceAll('/[^/]*$', "/") + child
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/_guessParamListFormat.nf'

def _guessParamListFormat(params) {
  if (!params.containsKey("param_list") || params.param_list == null) {
    "none"
  } else {
    def param_list = params.param_list

    if (param_list !instanceof String) {
      "asis"
    } else if (param_list.endsWith(".csv")) {
      "csv"
    } else if (param_list.endsWith(".json") || param_list.endsWith(".jsn")) {
      "json"
    } else if (param_list.endsWith(".yaml") || param_list.endsWith(".yml")) {
      "yaml"
    } else {
      "yaml_blob"
    }
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/IDChecker.nf'
class IDChecker {
  final def items = [] as Set

  @groovy.transform.WithWriteLock
  boolean observe(String item) {
    if (items.contains(item)) {
      return false
    } else {
      items << item
      return true
    }
  }

  @groovy.transform.WithReadLock
  boolean contains(String item) {
    return items.contains(item)
  }

  @groovy.transform.WithReadLock
  Set getItems() {
    return items.clone()
  }
}
// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/paramsToChannel.nf'
def paramsToChannel(params, config) {
  if (!viashChannelDeprecationWarningPrinted) {
    viashChannelDeprecationWarningPrinted = true
    System.err.println("Warning: paramsToChannel has deprecated in Viash 0.7.0. " +
                      "Please use a combination of channelFromParams and preprocessInputs.")
  }
  Channel.fromList(paramsToList(params, config))
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/paramsToList.nf'
viashChannelDeprecationWarningPrinted = false

def paramsToList(params, config) {
  if (!viashChannelDeprecationWarningPrinted) {
    viashChannelDeprecationWarningPrinted = true
    System.err.println("Warning: paramsToList has deprecated in Viash 0.7.0. " +
                      "Please use a combination of channelFromParams and preprocessInputs.")
  }
  // fetch default params from functionality
  def defaultArgs = config.functionality.allArguments
    .findAll { it.containsKey("default") }
    .collectEntries { [ it.plainName, it.default ] }

  // fetch overrides in params
  def paramArgs = config.functionality.allArguments
    .findAll { params.containsKey(it.plainName) }
    .collectEntries { [ it.plainName, params[it.plainName] ] }
  
  // check multi input params
  // objects should be closures and not functions, thanks to FunctionDef
  def multiParamFormat = _guessParamListFormat(params)

  def multiOptionFunctions = [ 
    "csv": {[it, readCsv(it)]},
    "json": {[it, readJson(it)]},
    "yaml": {[it, readYaml(it)]},
    "yaml_blob": {[null, readYamlBlob(it)]},
    "asis": {[null, it]},
    "none": {[null, [[:]]]}
  ]
  assert multiOptionFunctions.containsKey(multiParamFormat): 
    "Format of provided --param_list not recognised.\n" +
    "You can use '--param_list_format' to manually specify the format.\n" +
    "Found: '$multiParamFormat'. Expected: one of 'csv', 'json', 'yaml', 'yaml_blob', 'asis' or 'none'"

  // fetch multi param inputs
  def multiOptionFun = multiOptionFunctions.get(multiParamFormat)
  // todo: add try catch
  def multiOptionOut = multiOptionFun(params.containsKey("param_list") ? params.param_list : "")
  def paramList = multiOptionOut[1]
  def multiFile = multiOptionOut[0]

  // data checks
  assert paramList instanceof List: "--param_list should contain a list of maps"
  for (value in paramList) {
    assert value instanceof Map: "--param_list should contain a list of maps"
  }
  
  // combine parameters
  def processedParams = paramList.collect{ multiParam ->
    // combine params
    def combinedArgs = defaultArgs + paramArgs + multiParam

    if (workflow.stubRun) {
      // if stub run, explicitly add an id if missing
      combinedArgs = [id: "stub"] + combinedArgs
    } else {
      // else check whether required arguments exist
      config.functionality.allArguments
        .findAll { it.required }
        .forEach { par ->
          assert combinedArgs.containsKey(par.plainName): "Argument ${par.plainName} is required but does not have a value"
        }
    }
    
    // process arguments
    def inputs = config.functionality.allArguments
      .findAll{ par -> combinedArgs.containsKey(par.plainName) }
      .collectEntries { par ->
        // split on 'multiple_sep'
        if (par.multiple) {
          parData = combinedArgs[par.plainName]
          if (parData instanceof List) {
            parData = parData.collect{it instanceof String ? it.split(par.multiple_sep) : it }
          } else if (parData instanceof String) {
            parData = parData.split(par.multiple_sep)
          } else if (parData == null) {
            parData = []
          } else {
            parData = [ parData ]
          }
        } else {
          parData = [ combinedArgs[par.plainName] ]
        }

        // flatten
        parData = parData.flatten()

        // cast types
        if (par.type == "file" && ((par.direction != null ? par.direction : "input") == "input")) {
          parData = parData.collect{path ->
            if (path !instanceof String) {
              path
            } else if (multiFile) {
              file(_getChild(multiFile, path))
            } else {
              file(path)
            }
          }.flatten()
        } else if (par.type == "integer") {
          parData = parData.collect{it as Integer}
        } else if (par.type == "double") {
          parData = parData.collect{it as Double}
        } else if (par.type == "boolean" || par.type == "boolean_true" || par.type == "boolean_false") {
          parData = parData.collect{it as Boolean}
        }
        // simplify list to value if need be
        if (!par.multiple) {
          assert parData.size() == 1 : 
            "Error: argument ${par.plainName} has too many values.\n" +
            "  Expected amount: 1. Found: ${parData.size()}"
          parData = parData[0]
        }

        // return pair
        [ par.plainName, parData ]
      }
      // remove parameters which were explicitly set to null
      .findAll{ par -> par != null }
    }
    
  
  // check processed params
  processedParams.forEach { args ->
    assert args.containsKey("id"): "Each argument set should have an 'id'. Argument set: $args"
  }
  def ppIds = processedParams.collect{it.id}
  assert ppIds.size() == ppIds.unique().size() : "All argument sets should have unique ids. Detected ids: $ppIds"

  processedParams
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/preprocessInputs.nf'
// This helper file will be deprecated soon
preprocessInputsDeprecationWarningPrinted = false

def preprocessInputsDeprecationWarning() {
  if (!preprocessInputsDeprecationWarningPrinted) {
    preprocessInputsDeprecationWarningPrinted = true
    System.err.println("Warning: preprocessInputs() will be deprecated Viash 0.9.0.")
  }
}

/**
 * Process a list of Vdsl3 formatted parameters and apply a Viash config to them:
 *    - Gather default parameters from the Viash config and make 
 *      sure that they are correctly formatted (see applyConfig method).
 *    - Format the input parameters (also using the applyConfig method).
 *    - Apply the default parameter to the input parameters.
 *    - Do some assertions:
 *        ~ Check if the event IDs in the channel are unique.
 *  
 * @param params A list of parameter sets as Tuples. The first element of the tuples
 *                must be a unique id of the parameter set, and the second element 
 *                must contain the parameters themselves. Optional extra elements 
 *                of the tuples will be passed to the output as is.
 * @param config A Map of the Viash configuration. This Map can be generated from 
 *                the config file using the readConfig() function.           
 *
 * @return A list of processed parameters sets as tuples.
 */

private List<Tuple> _preprocessInputsList(List<Tuple> params, Map config) {
  // Get different parameter types (used throughout this function)
  def defaultArgs = config.functionality.allArguments
    .findAll { it.containsKey("default") }
    .collectEntries { [ it.plainName, it.default ] }

  // Apply config to default parameters
  def parsedDefaultValues = applyConfigToOneParameterSet(defaultArgs, config)

  // Apply config to input parameters
  def parsedInputParamSets = applyConfig(params, config)

  // Merge two parameter sets together
  def parsedArgs = parsedInputParamSets.collect({ parsedInputParamSet ->
    def id = parsedInputParamSet[0]
    def parValues = parsedInputParamSet[1]
    def passthrough = parsedInputParamSet.drop(2)
    def parValuesWithDefault = parsedDefaultValues + parValues
    [id, parValuesWithDefault] + passthrough
  })
  _checkUniqueIds(parsedArgs)

  return parsedArgs
}

/**
 * Generate a nextflow Workflow that allows processing a channel of 
 * Vdsl3 formatted events and apply a Viash config to them:
 *    - Gather default parameters from the Viash config and make 
 *      sure that they are correctly formatted (see applyConfig method).
 *    - Format the input parameters (also using the applyConfig method).
 *    - Apply the default parameter to the input parameters.
 *    - Do some assertions:
 *        ~ Check if the event IDs in the channel are unique.
 * 
 * The events in the channel are formatted as tuples, with the 
 * first element of the tuples being a unique id of the parameter set, 
 * and the second element containg the the parameters themselves.
 * Optional extra elements of the tuples will be passed to the output as is.
 *
 * @param args A map that must contain a 'config' key that points
 *              to a parsed config (see readConfig()). Optionally, a
 *              'key' key can be provided which can be used to create a unique
 *              name for the workflow process.
 *
 * @return A workflow that allows processing a channel of Vdsl3 formatted events
 * and apply a Viash config to them.
 */
def preprocessInputs(Map args) {
  preprocessInputsDeprecationWarning()
  
  wfKey = args.key != null ? args.key : "preprocessInputs"
  config = args.config
  workflow preprocessInputsInstance {
    take: 
    input_ch

    main:
    assert config instanceof Map : 
      "Error in preprocessInputs: config must be a map. " +
      "Expected class: Map. Found: config.getClass() is ${config.getClass()}"

    output_ch = input_ch
      | toSortedList
      | map { paramList -> _preprocessInputsList(paramList, config) }
      | flatMap
    emit:
    output_ch
  }

  return preprocessInputsInstance.cloneWithName(wfKey)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/runComponents.nf'
/**
 * Run a list of components on a stream of data.
 * 
 * @param components: list of Viash VDSL3 modules to run
 * @param fromState: a closure, a map or a list of keys to extract from the input data.
 *   If a closure, it will be called with the id, the data and the component config.
 * @param toState: a closure, a map or a list of keys to extract from the output data
 *   If a closure, it will be called with the id, the output data, the old state and the component config.
 * @param filter: filter function to apply to the input.
 *   It will be called with the id, the data and the component config.
 * @param id: id to use for the output data
 *   If a closure, it will be called with the id, the data and the component config.
 * @param auto: auto options to pass to the components
 *
 * @return: a workflow that runs the components
 **/
def runComponents(Map args) {
  assert args.components: "runComponents should be passed a list of components to run"

  def components_ = args.components
  if (components_ !instanceof List) {
    components_ = [ components_ ]
  }
  assert components_.size() > 0: "pass at least one component to runComponents"

  def fromState_ = args.fromState
  def toState_ = args.toState
  def filter_ = args.filter
  def id_ = args.id

  workflow runComponentsWf {
    take: input_ch
    main:

    // generate one channel per method
    out_chs = components_.collect{ comp_ ->
      def comp_config = comp_.config

      filter_ch = filter_
        ? input_ch | filter{tup ->
          filter_(tup[0], tup[1], comp_config)
        }
        : input_ch
      id_ch = id_
        ? filter_ch | map{tup ->
          // def new_id = id_(tup[0], tup[1], comp_config)
          def new_id = tup[0]
          if (id_ instanceof String) {
            new_id = id_
          } else if (id_ instanceof Closure) {
            new_id = id_(new_id, tup[1], comp_config)
          }
          [new_id] + tup.drop(1)
        }
        : filter_ch
      data_ch = id_ch | map{tup ->
          def new_data = tup[1]
          if (fromState_ instanceof Map) {
            new_data = fromState_.collectEntries{ key0, key1 ->
              [key0, new_data[key1]]
            }
          } else if (fromState_ instanceof List) {
            new_data = fromState_.collectEntries{ key ->
              [key, new_data[key]]
            }
          } else if (fromState_ instanceof Closure) {
            new_data = fromState_(tup[0], new_data, comp_config)
          }
          tup.take(1) + [new_data] + tup.drop(1)
        }
      out_ch = data_ch
        | comp_.run(
          auto: (args.auto ?: [:]) + [simplifyInput: false, simplifyOutput: false]
        )
      post_ch = toState_
        ? out_ch | map{tup ->
          def output = tup[1]
          def old_state = tup[2]
          if (toState_ instanceof Map) {
            new_state = old_state + toState_.collectEntries{ key0, key1 ->
              [key0, output[key1]]
            }
          } else if (toState_ instanceof List) {
            new_state = old_state + toState_.collectEntries{ key ->
              [key, output[key]]
            }
          } else if (toState_ instanceof Closure) {
            new_state = toState_(tup[0], output, old_state, comp_config)
          }
          [tup[0], new_state] + tup.drop(3)
        }
        : out_ch
      
      post_ch
    }

    // mix all results
    output_ch =
      (out_chs.size == 1)
        ? out_chs[0]
        : out_chs[0].mix(*out_chs.drop(1))

    emit: output_ch
  }

  return runComponentsWf
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/safeJoin.nf'
def safeJoin(targetChannel, sourceChannel, key) {
  def sourceIDs = new IDChecker()

  def sourceCheck = sourceChannel
    | map { tup ->
      sourceIDs.observe(tup[0])
      tup
    }
  def targetCheck = targetChannel
    | map { tup ->
      def id = tup[0]
      // def id = tup[1].containsKey("_meta").containsKey("join_id") ? tup[1]._meta.join_id : tup[0]
      
      if (!sourceIDs.contains(id)) {
        error (
          "Error in module '${key}' when merging output with original state.\n" +
          "  Reason: output with id '${id}' could not be joined with source channel.\n" +
          "    If the IDs in the output channel differ from the input channel,\n" + 
          "    please set `tup[1]._meta.join_id to the original ID.\n" +
          "  Original IDs in input channel: ['${sourceIDs.getItems().join("', '")}'].\n" + 
          "  Unexpected ID in the output channel: '${id}.\n" +
          "  Example input event: [\"id\", [input: file(...)]],\n" +
          "  Example output event: [\"newid\", [output: file(...), _meta: [join_id: \"id\"]]]"
        )
      }

      tup
    }
  targetCheck.join(sourceCheck)
}
// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/_splitParams.nf'
/**
 * Split parameters for arguments that accept multiple values using their separator
 *
 * @param paramList A Map containing parameters to split.
 * @param config A Map of the Viash configuration. This Map can be generated from the config file
 *               using the readConfig() function.
 *
 * @return A Map of parameters where the parameter values have been split into a list using
 *         their seperator.
 */
Map<String, Object> _splitParams(Map<String, Object> parValues, Map config){
  def parsedParamValues = parValues.collectEntries { parName, parValue ->
    def parameterSettings = config.functionality.allArguments.find({it.plainName == parName})

    if (!parameterSettings) {
      // if argument is not found, do not alter 
      return [parName, parValue]
    }
    if (parameterSettings.multiple) { // Check if parameter can accept multiple values
      if (parValue instanceof Collection) {
          parValue = parValue.collect{it instanceof String ? it.split(parameterSettings.multiple_sep) : it }
      } else if (parValue instanceof String) {
          parValue = parValue.split(parameterSettings.multiple_sep)
      } else if (parValue == null) {
          parValue = []
      } else {
          parValue = [ parValue ]
      }
      parValue = parValue.flatten()
    }
    // For all parameters check if multiple values are only passed for
    // arguments that allow it. Quietly simplify lists of length 1.
    if (!parameterSettings.multiple && parValue instanceof Collection) {
      assert parValue.size() == 1 : 
      "Error: argument ${parName} has too many values.\n" +
      "  Expected amount: 1. Found: ${parValue.size()}"
      parValue = parValue[0]
    }
    [parName, parValue]
  }
  return parsedParamValues
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/channel/viashChannel.nf'

def viashChannel(params, config) {
  if (!viashChannelDeprecationWarningPrinted) {
    viashChannelDeprecationWarningPrinted = true
    System.err.println("Warning: viashChannel has deprecated in Viash 0.7.0. " +
                      "Please use a combination of channelFromParams and preprocessInputs.")
  }
  paramsToChannel(params, config)
    | map{tup -> [tup.id, tup]}
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/config/addGlobalParams.nf'
def addGlobalArguments(config) {
  def localConfig = [
    "functionality" : [
      "argument_groups": [
        [
          "name": "Nextflow input-output arguments",
          "description": "Input/output parameters for Nextflow itself. Please note that both publishDir and publish_dir are supported but at least one has to be configured.",
          "arguments" : [
            [
              'name': '--publish_dir',
              'required': true,
              'type': 'string',
              'description': 'Path to an output directory.',
              'example': 'output/',
              'multiple': false
            ],
            [
              'name': '--param_list',
              'required': false,
              'type': 'string',
              'description': '''Allows inputting multiple parameter sets to initialise a Nextflow channel. A `param_list` can either be a list of maps, a csv file, a json file, a yaml file, or simply a yaml blob.
              |
              |* A list of maps (as-is) where the keys of each map corresponds to the arguments of the pipeline. Example: in a `nextflow.config` file: `param_list: [ ['id': 'foo', 'input': 'foo.txt'], ['id': 'bar', 'input': 'bar.txt'] ]`.
              |* A csv file should have column names which correspond to the different arguments of this pipeline. Example: `--param_list data.csv` with columns `id,input`.
              |* A json or a yaml file should be a list of maps, each of which has keys corresponding to the arguments of the pipeline. Example: `--param_list data.json` with contents `[ {'id': 'foo', 'input': 'foo.txt'}, {'id': 'bar', 'input': 'bar.txt'} ]`.
              |* A yaml blob can also be passed directly as a string. Example: `--param_list "[ {'id': 'foo', 'input': 'foo.txt'}, {'id': 'bar', 'input': 'bar.txt'} ]"`.
              |
              |When passing a csv, json or yaml file, relative path names are relativized to the location of the parameter file. No relativation is performed when `param_list` is a list of maps (as-is) or a yaml blob.'''.stripMargin(),
              'example': 'my_params.yaml',
              'multiple': false,
              'hidden': true
            ],
          ]
        ]
      ]
    ]
  ]

  return processConfig(_mergeMap(config, localConfig))
}

def _mergeMap(Map lhs, Map rhs) {
  return rhs.inject(lhs.clone()) { map, entry ->
    if (map[entry.key] instanceof Map && entry.value instanceof Map) {
      map[entry.key] = _mergeMap(map[entry.key], entry.value)
    } else if (map[entry.key] instanceof Collection && entry.value instanceof Collection) {
      map[entry.key] += entry.value
    } else {
      map[entry.key] = entry.value
    }
    return map
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/config/generateHelp.nf'
def _generateArgumentHelp(param) {
  // alternatives are not supported
  // def names = param.alternatives ::: List(param.name)

  def unnamedProps = [
    ["required parameter", param.required],
    ["multiple values allowed", param.multiple],
    ["output", param.direction.toLowerCase() == "output"],
    ["file must exist", param.type == "file" && param.must_exist]
  ].findAll{it[1]}.collect{it[0]}
  
  def dflt = null
  if (param.default != null) {
    if (param.default instanceof List) {
      dflt = param.default.join(param.multiple_sep != null ? param.multiple_sep : ", ")
    } else {
      dflt = param.default.toString()
    }
  }
  def example = null
  if (param.example != null) {
    if (param.example instanceof List) {
      example = param.example.join(param.multiple_sep != null ? param.multiple_sep : ", ")
    } else {
      example = param.example.toString()
    }
  }
  def min = param.min?.toString()
  def max = param.max?.toString()

  def escapeChoice = { choice ->
    def s1 = choice.replaceAll("\\n", "\\\\n")
    def s2 = s1.replaceAll("\"", """\\\"""")
    s2.contains(",") || s2 != choice ? "\"" + s2 + "\"" : s2
  }
  def choices = param.choices == null ? 
    null : 
    "[ " + param.choices.collect{escapeChoice(it.toString())}.join(", ") + " ]"

  def namedPropsStr = [
    ["type", ([param.type] + unnamedProps).join(", ")],
    ["default", dflt],
    ["example", example],
    ["choices", choices],
    ["min", min],
    ["max", max]
  ]
    .findAll{it[1]}
    .collect{"\n        " + it[0] + ": " + it[1].replaceAll("\n", "\\n")}
    .join("")
  
  def descStr = param.description == null ?
    "" :
    _paragraphWrap("\n" + param.description.trim(), 80 - 8).join("\n        ")
  
  "\n    --" + param.plainName +
    namedPropsStr +
    descStr
}

// Based on Helper.generateHelp() in Helper.scala
def _generateHelp(config) {
  def fun = config.functionality

  // PART 1: NAME AND VERSION
  def nameStr = fun.name + 
    (fun.version == null ? "" : " " + fun.version)

  // PART 2: DESCRIPTION
  def descrStr = fun.description == null ? 
    "" :
    "\n\n" + _paragraphWrap(fun.description.trim(), 80).join("\n")

  // PART 3: Usage
  def usageStr = fun.usage == null ? 
    "" :
    "\n\nUsage:\n" + fun.usage.trim()

  // PART 4: Options
  def argGroupStrs = fun.allArgumentGroups.collect{argGroup ->
    def name = argGroup.name
    def descriptionStr = argGroup.description == null ?
      "" :
      "\n    " + _paragraphWrap(argGroup.description.trim(), 80-4).join("\n    ") + "\n"
    def arguments = argGroup.arguments.collect{arg -> 
      arg instanceof String ? fun.allArguments.find{it.plainName == arg} : arg
    }.findAll{it != null}
    def argumentStrs = arguments.collect{param -> _generateArgumentHelp(param)}
    
    "\n\n$name:" +
      descriptionStr +
      argumentStrs.join("\n")
  }

  // FINAL: combine
  def out = nameStr + 
    descrStr +
    usageStr + 
    argGroupStrs.join("")

  return out
}

// based on Format._paragraphWrap
def _paragraphWrap(str, maxLength) {
  def outLines = []
  str.split("\n").each{par ->
    def words = par.split("\\s").toList()

    def word = null
    def line = words.pop()
    while(!words.isEmpty()) {
      word = words.pop()
      if (line.length() + word.length() + 1 <= maxLength) {
        line = line + " " + word
      } else {
        outLines.add(line)
        line = word
      }
    }
    if (words.isEmpty()) {
      outLines.add(line)
    }
  }
  return outLines
}

def helpMessage(config) {
  if (params.containsKey("help") && params.help) {
    def mergedConfig = addGlobalArguments(config)
    def helpStr = _generateHelp(mergedConfig)
    println(helpStr)
    exit 0
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/config/processConfig.nf'
def processConfig(config) {
  // TODO: assert .functionality etc.
  if (config.functionality.inputs) {
    System.err.println("Warning: .functionality.inputs is deprecated. Please use .functionality.arguments instead.")
  }
  if (config.functionality.outputs) {
    System.err.println("Warning: .functionality.outputs is deprecated. Please use .functionality.arguments instead.")
  }

  // set defaults for inputs
  config.functionality.inputs = 
    (config.functionality.inputs != null ? config.functionality.inputs : []).collect{arg ->
      arg.type = arg.type != null ? arg.type : "file"
      arg.direction = "input"
      _processArgument(arg)
    }
  // set defaults for outputs
  config.functionality.outputs = 
    (config.functionality.outputs != null ? config.functionality.outputs : []).collect{arg ->
      arg.type = arg.type != null ? arg.type : "file"
      arg.direction = "output"
      _processArgument(arg)
    }
  // set defaults for arguments
  config.functionality.arguments = 
    (config.functionality.arguments != null ? config.functionality.arguments : []).collect{arg ->
      _processArgument(arg)
    }
  // set defaults for argument_group arguments
  config.functionality.argument_groups =
    (config.functionality.argument_groups != null ? config.functionality.argument_groups : []).collect{grp ->
      grp.arguments = (grp.arguments != null ? grp.arguments : []).collect{arg ->
        arg instanceof String ? arg.replaceAll("^-*", "") : _processArgument(arg)
      }
      grp
    }

  // create combined arguments list
  config.functionality.allArguments = 
    config.functionality.inputs +
    config.functionality.outputs +
    config.functionality.arguments +
    config.functionality.argument_groups.collectMany{ group ->
      group.arguments.findAll{ it !instanceof String }
    }
  
  // add missing argument groups (based on Functionality::allArgumentGroups())
  def argGroups = config.functionality.argument_groups
  def inputGroup = _processArgumentGroup(argGroups, "Inputs", config.functionality.inputs)
  def outputGroup = _processArgumentGroup(argGroups, "Outputs", config.functionality.outputs)
  def defaultGroup = _processArgumentGroup(argGroups, "Arguments", config.functionality.arguments)
  def groupsFiltered = argGroups.findAll(gr -> !(["Inputs", "Outputs", "Arguments"].contains(gr.name)))
  config.functionality.allArgumentGroups = inputGroup + outputGroup + defaultGroup + groupsFiltered

  config
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/config/readConfig.nf'

def readConfig(file) {
  def config = readYaml(file != null ? file : "$projectDir/config.vsh.yaml")
  processConfig(config)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/functions/collectTraces.nf'
class CustomTraceObserver implements nextflow.trace.TraceObserver {
  List traces

  CustomTraceObserver(List traces) {
    this.traces = traces
  }

  @Override
  void onProcessComplete(nextflow.processor.TaskHandler handler, nextflow.trace.TraceRecord trace) {
    def trace2 = trace.store.clone()
    trace2.script = null
    traces.add(trace2)
  }

  @Override
  void onProcessCached(nextflow.processor.TaskHandler handler, nextflow.trace.TraceRecord trace) {
    def trace2 = trace.store.clone()
    trace2.script = null
    traces.add(trace2)
  }
}

def collectTraces() {
  def traces = Collections.synchronizedList([])

  // add custom trace observer which stores traces in the traces object
  session.observers.add(new CustomTraceObserver(traces))

  traces
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/functions/getPublishDir.nf'
def getPublishDir() {
  return params.containsKey("publish_dir") ? params.publish_dir : 
    params.containsKey("publishDir") ? params.publishDir : 
    null
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/functions/getRootDir.nf'

// Recurse upwards until we find a '.build.yaml' file
def _findBuildYamlFile(path) {
  def child = path.resolve(".build.yaml")
  if (java.nio.file.Files.isDirectory(path) && java.nio.file.Files.exists(child)) {
    return child
  } else {
    def parent = path.getParent()
    if (parent == null) {
      return null
    } else {
      return _findBuildYamlFile(parent)
    }
  }
}

// get the root of the target folder
def getRootDir() {
  def dir = _findBuildYamlFile(projectDir.toAbsolutePath())
  assert dir != null: "Could not find .build.yaml in the folder structure"
  dir.getParent()
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/functions/iterateMap.nf'
/**
  * Recursively apply a function over the leaves of an object.
  * @param obj The object to iterate over.
  * @param fun The function to apply to each value.
  * @return The object with the function applied to each value.
  */
def iterateMap(obj, fun) {
  if (obj instanceof List && obj !instanceof String) {
    return obj.collect{item ->
      iterateMap(item, fun)
    }
  } else if (obj instanceof Map) {
    return obj.collectEntries{key, item ->
      [key.toString(), iterateMap(item, fun)]
    }
  } else {
    return fun(obj)
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/functions/niceView.nf'
/**
  * A view for printing the event of each channel as a YAML blob.
  * This is useful for debugging.
  */
def niceView() {
  workflow niceViewWf {
    take: input
    main:
      output = input
        | view{toYamlBlob(it)}
    emit: output
  }
  return niceViewWf
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readCsv.nf'

def readCsv(file_path) {
  def output = []
  def inputFile = file_path !instanceof Path ? file(file_path) : file_path

  // todo: allow escaped quotes in string
  // todo: allow single quotes?
  def splitRegex = java.util.regex.Pattern.compile(''',(?=(?:[^"]*"[^"]*")*[^"]*$)''')
  def removeQuote = java.util.regex.Pattern.compile('''"(.*)"''')

  def br = java.nio.file.Files.newBufferedReader(inputFile)

  def row = -1
  def header = null
  while (br.ready() && header == null) {
    def line = br.readLine()
    row++
    if (!line.startsWith("#")) {
      header = splitRegex.split(line, -1).collect{field ->
        m = removeQuote.matcher(field)
        m.find() ? m.replaceFirst('$1') : field
      }
    }
  }
  assert header != null: "CSV file should contain a header"

  while (br.ready()) {
    def line = br.readLine()
    row++
    if (line == null) {
      br.close()
      break
    }

    if (!line.startsWith("#")) {
      def predata = splitRegex.split(line, -1)
      def data = predata.collect{field ->
        if (field == "") {
          return null
        }
        m = removeQuote.matcher(field)
        if (m.find()) {
          return m.replaceFirst('$1')
        } else {
          return field
        }
      }
      assert header.size() == data.size(): "Row $row should contain the same number as fields as the header"
      
      def dataMap = [header, data].transpose().collectEntries().findAll{it.value != null}
      output.add(dataMap)
    }
  }

  output
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readJsonBlob.nf'
def readJsonBlob(str) {
  def jsonSlurper = new groovy.json.JsonSlurper()
  jsonSlurper.parseText(str)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readJson.nf'
def readJson(file_path) {
  def inputFile = file_path !instanceof Path ? file(file_path) : file_path
  def jsonSlurper = new groovy.json.JsonSlurper()
  jsonSlurper.parse(inputFile)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readTaggedYaml.nf'
// Custom constructor to modify how certain objects are parsed from YAML
class CustomConstructor extends org.yaml.snakeyaml.constructor.Constructor {
  Path root

  class ConstructPath extends org.yaml.snakeyaml.constructor.AbstractConstruct {
    public Object construct(org.yaml.snakeyaml.nodes.Node node) {
      String filename = (String) constructScalar(node);
      if (root != null) {
        return root.resolve(filename);
      }
      return java.nio.file.Paths.get(filename);
    }
  }

  CustomConstructor(org.yaml.snakeyaml.LoaderOptions options, Path root) {
    super(options)
    this.root = root
    // Handling !file tag and parse it back to a File type
    this.yamlConstructors.put(new org.yaml.snakeyaml.nodes.Tag("!file"), new ConstructPath())
  }
}

def readTaggedYaml(Path path) {
  def options = new org.yaml.snakeyaml.LoaderOptions()
  def constructor = new CustomConstructor(options, path.getParent())
  def yaml = new org.yaml.snakeyaml.Yaml(constructor)
  return yaml.load(path.text)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readYamlBlob.nf'
def readYamlBlob(str) {
  def yamlSlurper = new org.yaml.snakeyaml.Yaml()
  yamlSlurper.load(str)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/readYaml.nf'
def readYaml(file_path) {
  def inputFile = file_path !instanceof Path ? file(file_path) : file_path
  def yamlSlurper = new org.yaml.snakeyaml.Yaml()
  yamlSlurper.load(inputFile)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/toJsonBlob.nf'
String toJsonBlob(Map data) {
  return groovy.json.JsonOutput.toJson(data)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/toTaggedYamlBlob.nf'
// Custom representer to modify how certain objects are represented in YAML
class CustomRepresenter extends org.yaml.snakeyaml.representer.Representer {
  class RepresentPath implements org.yaml.snakeyaml.representer.Represent {
    public org.yaml.snakeyaml.nodes.Node representData(Object data) {
      Path file = (Path) data;
      String value = file.getFileName();
      def tag = new org.yaml.snakeyaml.nodes.Tag("!file");
      return representScalar(tag, value);
    }
  }
  CustomRepresenter(org.yaml.snakeyaml.DumperOptions options) {
    super(options)
    this.representers.put(Path, new RepresentPath())
  }
}

String toTaggedYamlBlob(Map data) {
  def options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  def representer = new CustomRepresenter(options)
  def yaml = new org.yaml.snakeyaml.Yaml(representer, options)
  return yaml.dump(data)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/toYamlBlob.nf'
String toYamlBlob(Map data) {
  def options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  options.setPrettyFlow(true)
  def yaml = new org.yaml.snakeyaml.Yaml(options)
  def cleanData = iterateMap(data, {it.toString()})
  return yaml.dump(data)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/writeJson.nf'
void writeJson(data, file) {
  assert data: "writeJson: data should not be null"
  assert file: "writeJson: file should not be null"
  file.write(toJsonBlob(data))
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/readwrite/writeYaml.nf'
void writeYaml(data, file) {
  assert data: "writeYaml: data should not be null"
  assert file: "writeYaml: file should not be null"
  file.write(toYamlBlob(data))
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/states/findStates.nf'
def findStates(Map params, Map config) {
  // TODO: do a deep clone of config
  def auto_config = config.clone()
  auto_config.functionality = auto_config.functionality.clone()
  // override arguments
  auto_config.functionality.argument_groups = []
  auto_config.functionality.arguments = [
    [
      type: "file",
      name: "--input_states",
      example: "/path/to/input/directory/**/state.yaml",
      description: "Path to input directory containing the datasets to be integrated.",
      required: true,
      multiple: true,
      multiple_sep: ";"
    ],
    [
      type: "string",
      name: "--filter",
      example: "foo/.*/state.yaml",
      description: "Regex to filter state files by path.",
      required: false
    ],
    // to do: make this a yaml blob?
    [
      type: "string",
      name: "--rename_keys",
      example: ["newKey1:oldKey1", "newKey2:oldKey2"],
      description: "Rename keys in the detected input files. This is useful if the input files do not match the set of input arguments of the workflow.",
      required: false,
      multiple: true,
      multiple_sep: ","
    ],
    [
      type: "string",
      name: "--settings",
      example: '{"output_dataset": "dataset.h5ad", "k": 10}',
      description: "Global arguments as a JSON glob to be passed to all components.",
      required: false
    ]
  ]

  // run auto config through processConfig once more
  auto_config = processConfig(auto_config)

  workflow findStatesWf {
    helpMessage(auto_config)

    output_ch = 
      channelFromParams(params, auto_config)
        | flatMap { autoId, args ->

          def globalSettings = args.settings ? readYamlBlob(args.settings) : [:]

          // look for state files in input dir
          def stateFiles = args.input_states

          // filter state files by regex
          if (args.filter) {
            stateFiles = stateFiles.findAll{ stateFile ->
              def stateFileStr = stateFile.toString()
              def matcher = stateFileStr =~ args.filter
              matcher.matches()}
          }

          // read in states
          def states = stateFiles.collect { stateFile ->
            def state_ = readTaggedYaml(stateFile)
            [state_.id, state_]
          }

          // construct renameMap
          if (args.rename_keys) {
            def renameMap = args.rename_keys.collectEntries{renameString ->
              def split = renameString.split(":")
              assert split.size() == 2: "Argument 'rename' should be of the form 'newKey:oldKey,newKey:oldKey'"
              split
            }

            // rename keys in state, only let states through which have all keys
            // also add global settings
            states = states.collectMany{id, state ->
              def newState = [:]

              for (key in renameMap.keySet()) {
                def origKey = renameMap[key]
                if (!(state.containsKey(origKey))) {
                  return []
                }
                newState[key] = state[origKey]
              }

              [[id, globalSettings + newState]]
            }
          }

          states
        }
    emit:
    output_ch
  }

  return findStatesWf
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/states/publishStates.nf'
def collectFiles(obj) {
  if (obj instanceof java.io.File || obj instanceof Path)  {
    return [obj]
  } else if (obj instanceof List && obj !instanceof String) {
    return obj.collectMany{item ->
      collectFiles(item)
    }
  } else if (obj instanceof Map) {
    return obj.collectMany{key, item ->
      collectFiles(item)
    }
  } else {
    return []
  }
}

/**
 * Recurse through a state and collect all input files and their target output filenames.
 * @param obj The state to recurse through.
 * @param prefix The prefix to prepend to the output filenames.
 */
def collectInputOutputPaths(obj, prefix) {
  if (obj instanceof File || obj instanceof Path)  {
    def path = obj instanceof Path ? obj : obj.toPath()
    def ext = path.getFileName().find("\\.[^\\.]+\$") ?: ""
    def newFilename = prefix + ext
    return [[obj, newFilename]]
  } else if (obj instanceof List && obj !instanceof String) {
    return obj.withIndex().collectMany{item, ix ->
      collectInputOutputPaths(item, prefix + "_" + ix)
    }
  } else if (obj instanceof Map) {
    return obj.collectMany{key, item ->
      collectInputOutputPaths(item, prefix + "." + key)
    }
  } else {
    return []
  }
}

def publishStates(Map args) {
  def key_ = args.get("key")
  
  workflow publishStatesWf {
    take: input_ch
    main:
      input_ch
        | map { tup ->
          def id_ = tup[0]
          def state_ = tup[1]

          // the input files and the target output filenames
          def inputOutputFiles_ = collectInputOutputPaths(state_, id_ + "." + key_).transpose()
          def inputFiles_ = inputOutputFiles_[0]
          def outputFiles_ = inputOutputFiles_[1]

          // convert state to yaml blob
          def yamlBlob_ = toTaggedYamlBlob([id: id_] + state_)

          // adds a leading dot to the id (after any folder names)
          // example: foo -> .foo, foo/bar -> foo/.bar
          def idWithDot_ = id_.replaceAll("^(.+/)?([^/]+)", "\$1.\$2")
          def yamlFile = '$id.$key.state.yaml'
            .replaceAll('\\$id', idWithDot_)
            .replaceAll('\\$key', key_)

          [id_, yamlBlob_, yamlFile, inputFiles_, outputFiles_]
        }
        | publishStatesProc
    emit: input_ch
  }
  return publishStatesWf
}
process publishStatesProc {
  // todo: check publishpath?
  publishDir path: "${getPublishDir()}/", mode: "copy"
  tag "$id"
  input:
    tuple val(id), val(yamlBlob), val(yamlFile), path(inputFiles, stageAs: "_inputfile?/*"), val(outputFiles)
  output:
    tuple val(id), path{[yamlFile] + outputFiles}
  script:
  def copyCommands = [
    inputFiles instanceof List ? inputFiles : [inputFiles],
    outputFiles instanceof List ? outputFiles : [outputFiles]
  ]
    .transpose()
    .collectMany{infile, outfile ->
      if (infile.toString() != outfile.toString()) {
        ["cp -r '${infile.toString()}' '${outfile.toString()}'"]
      } else {
        // no need to copy if infile is the same as outfile
        []
      }
    }
  """
  mkdir -p "\$(dirname '${yamlFile}')"
  echo "Storing state as yaml"
  echo '${yamlBlob}' > '${yamlFile}'
  echo "Copying output files to destination folder"
  ${copyCommands.join("\n  ")}
  """
}


// this assumes that the state contains no other values other than those specified in the config
def publishStatesByConfig(Map args) {
  def key_ = args.get("key")
  def config = args.get("config")
  
  workflow publishStatesSimpleWf {
    take: input_ch
    main:
      input_ch
        | map { tup ->
          def id_ = tup[0]
          def state_ = tup[1] // e.g. [output: new File("myoutput.h5ad"), k: 10]
          def origState_ = tup[2] // e.g. [output: '$id.$key.foo.h5ad']

          // the processed state is a list of [key, value, srcPath, destPath] tuples, where
          //   - key, value is part of the state to be saved to disk
          //   - srcPath and destPath are lists of files to be copied from src to dest
          def processedState =
            config.functionality.allArguments
              .findAll { it.direction == "output" }
              .collectMany { par ->
                def plainName_ = par.plainName
                // if the state does not contain the key, it's an
                // optional argument for which the component did 
                // not generate any output
                if (!state_.containsKey(plainName_)) {
                  return []
                }
                def value = state_[plainName_]
                // if the parameter is not a file, it should be stored
                // in the state as-is, but is not something that needs 
                // to be copied from the source path to the dest path
                if (par.type != "file") {
                  return [[key: plainName_, value: value, srcPath: [], destPath: []]]
                }
                // if the orig state does not contain this filename,
                // it's an optional argument for which the user specified
                // that it should not be returned as a state
                if (!origState_.containsKey(plainName_)) {
                  return []
                }
                def filenameTemplate = origState_[plainName_]
                // if the pararameter is multiple: true, fetch the template
                if (par.multiple && filenameTemplate instanceof List) {
                  filenameTemplate = filenameTemplate[0]
                }
                // instantiate the template
                filename = filenameTemplate
                  .replaceAll('\\$id', id_)
                  .replaceAll('\\$key', key_)
                if (par.multiple) {
                  // if the parameter is multiple: true, the filename
                  // should contain a wildcard '*' that is replaced with
                  // the index of the file
                  assert filename.contains("*") : "Module '${key_}' id '${id_}': Multiple output files specified, but no wildcard '*' in the filename: ${filename}"
                  def outputPerFile = value.withIndex().collect{ val, ix ->
                    def destPath = filename.replace("*", ix.toString())
                    def destFile = java.nio.file.Paths.get(destPath)
                    def srcPath = val instanceof File ? val.toPath() : val
                    [value: destFile, srcPath: srcPath, destPath: destPath]
                  }
                  def transposedOutputs = ["value", "srcPath", "destPath"].collectEntries{ key -> 
                    [key, outputPerFile.collect{dic -> dic[key]}]
                  }
                  return [[key: plainName_] + transposedOutputs]
                } else {
                  def destFile = java.nio.file.Paths.get(filename)
                  def srcPath = value instanceof File ? value.toPath() : value
                  return [[key: plainName_, value: destFile, srcPath: [srcPath], destPath: [filename]]]
                }
              }
          
          def updatedState_ = processedState.collectEntries{[it.key, it.value]}
          def inputFiles_ = processedState.collectMany{it.srcPath}
          def outputFiles_ = processedState.collectMany{it.destPath}
          
          // convert state to yaml blob
          def yamlBlob_ = toTaggedYamlBlob([id: id_] + updatedState_)

          // adds a leading dot to the id (after any folder names)
          // example: foo -> .foo, foo/bar -> foo/.bar
          // TODO: allow defining the state.yaml template
          def idWithDot_ = id_.replaceAll("^(.+/)?([^/]+)", "\$1.\$2")
          def yamlFile = '$id.$key.state.yaml'
            .replaceAll('\\$id', idWithDot_)
            .replaceAll('\\$key', key_)

          [id_, yamlBlob_, yamlFile, inputFiles_, outputFiles_]
        }
        | publishStatesProc
    emit: input_ch
  }
  return publishStatesSimpleWf
}
// helper file: 'src/main/resources/io/viash/platforms/nextflow/states/setState.nf'
def setState(fun) {
  assert fun instanceof Closure || fun instanceof Map || fun instanceof List :
    "Error in setState: Expected process argument to be a Closure, a Map, or a List. Found: class ${fun.getClass()}"

  // if fun is a List, convert to map
  if (fun instanceof List) {
    // check whether fun is a list[string]
    assert fun.every{it instanceof CharSequence} : "Error in setState: argument is a List, but not all elements are Strings"
    fun = fun.collectEntries{[it, it]}
  }

  // if fun is a map, convert to closure
  if (fun instanceof Map) {
    // check whether fun is a map[string, string]
    assert fun.values().every{it instanceof CharSequence} : "Error in setState: argument is a Map, but not all values are Strings"
    assert fun.keySet().every{it instanceof CharSequence} : "Error in setState: argument is a Map, but not all keys are Strings"
    def funMap = fun.clone()
    // turn the map into a closure to be used later on
    fun = { id_, state_ ->
      assert state_ instanceof Map : "Error in setState: the state is not a Map"
      funMap.collectMany{newkey, origkey ->
        if (state_.containsKey(origkey)) {
          [[newkey, state_[origkey]]]
        } else {
          []
        }
      }.collectEntries()
    }
  }

  map { tup ->
    def id = tup[0]
    def state = tup[1]
    def unfilteredState = fun(id, state)
    def newState = unfilteredState.findAll{key, val -> val != null}
    [id, newState] + tup.drop(2)
  }
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/processAuto.nf'
// TODO: unit test processAuto
def processAuto(Map auto) {
  // remove null values
  auto = auto.findAll{k, v -> v != null}

  // check for unexpected keys
  def expectedKeys = ["simplifyInput", "simplifyOutput", "transcript", "publish"]
  def unexpectedKeys = auto.keySet() - expectedKeys
  assert unexpectedKeys.isEmpty(), "unexpected keys in auto: '${unexpectedKeys.join("', '")}'"

  // check auto.simplifyInput
  assert auto.simplifyInput instanceof Boolean, "auto.simplifyInput must be a boolean"

  // check auto.simplifyOutput
  assert auto.simplifyOutput instanceof Boolean, "auto.simplifyOutput must be a boolean"

  // check auto.transcript
  assert auto.transcript instanceof Boolean, "auto.transcript must be a boolean"

  // check auto.publish
  assert auto.publish instanceof Boolean || auto.publish == "state", "auto.publish must be a boolean or 'state'"

  return auto.subMap(expectedKeys)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/processDirectives.nf'
def assertMapKeys(map, expectedKeys, requiredKeys, mapName) {
  assert map instanceof Map : "Expected argument '$mapName' to be a Map. Found: class ${map.getClass()}"
  map.forEach { key, val -> 
    assert key in expectedKeys : "Unexpected key '$key' in ${mapName ? mapName + " " : ""}map"
  }
  requiredKeys.forEach { requiredKey -> 
    assert map.containsKey(requiredKey) : "Missing required key '$key' in ${mapName ? mapName + " " : ""}map"
  }
}

// TODO: unit test processDirectives
def processDirectives(Map drctv) {
  // remove null values
  drctv = drctv.findAll{k, v -> v != null}

  // check for unexpected keys
  def expectedKeys = [
    "accelerator", "afterScript", "beforeScript", "cache", "conda", "container", "containerOptions", "cpus", "disk", "echo", "errorStrategy", "executor", "machineType", "maxErrors", "maxForks", "maxRetries", "memory", "module", "penv", "pod", "publishDir", "queue", "label", "scratch", "storeDir", "stageInMode", "stageOutMode", "tag", "time"
  ]
  def unexpectedKeys = drctv.keySet() - expectedKeys
  assert unexpectedKeys.isEmpty() : "Unexpected keys in process directive: '${unexpectedKeys.join("', '")}'"

  /* DIRECTIVE accelerator
    accepted examples:
    - [ limit: 4, type: "nvidia-tesla-k80" ]
  */
  if (drctv.containsKey("accelerator")) {
    assertMapKeys(drctv["accelerator"], ["type", "limit", "request", "runtime"], [], "accelerator")
  }

  /* DIRECTIVE afterScript
    accepted examples:
    - "source /cluster/bin/cleanup"
  */
  if (drctv.containsKey("afterScript")) {
    assert drctv["afterScript"] instanceof CharSequence
  }

  /* DIRECTIVE beforeScript
    accepted examples:
    - "source /cluster/bin/setup"
  */
  if (drctv.containsKey("beforeScript")) {
    assert drctv["beforeScript"] instanceof CharSequence
  }

  /* DIRECTIVE cache
    accepted examples:
    - true
    - false
    - "deep"
    - "lenient"
  */
  if (drctv.containsKey("cache")) {
    assert drctv["cache"] instanceof CharSequence || drctv["cache"] instanceof Boolean
    if (drctv["cache"] instanceof CharSequence) {
      assert drctv["cache"] in ["deep", "lenient"] : "Unexpected value for cache"
    }
  }

  /* DIRECTIVE conda
    accepted examples:
    - "bwa=0.7.15"
    - "bwa=0.7.15 fastqc=0.11.5"
    - ["bwa=0.7.15", "fastqc=0.11.5"]
  */
  if (drctv.containsKey("conda")) {
    if (drctv["conda"] instanceof List) {
      drctv["conda"] = drctv["conda"].join(" ")
    }
    assert drctv["conda"] instanceof CharSequence
  }

  /* DIRECTIVE container
    accepted examples:
    - "foo/bar:tag"
    - [ registry: "reg", image: "im", tag: "ta" ]
      is transformed to "reg/im:ta"
    - [ image: "im" ] 
      is transformed to "im:latest"
  */
  if (drctv.containsKey("container")) {
    assert drctv["container"] instanceof Map || drctv["container"] instanceof CharSequence
    if (drctv["container"] instanceof Map) {
      def m = drctv["container"]
      assertMapKeys(m, [ "registry", "image", "tag" ], ["image"], "container")
      def part1 = 
        System.getenv('OVERRIDE_CONTAINER_REGISTRY') ? System.getenv('OVERRIDE_CONTAINER_REGISTRY') + "/" : 
        params.containsKey("override_container_registry") ? params["override_container_registry"] + "/" : // todo: remove?
        m.registry ? m.registry + "/" : 
        ""
      def part2 = m.image
      def part3 = m.tag ? ":" + m.tag : ":latest"
      drctv["container"] = part1 + part2 + part3
    }
  }

  /* DIRECTIVE containerOptions
    accepted examples:
    - "--foo bar"
    - ["--foo bar", "-f b"]
  */
  if (drctv.containsKey("containerOptions")) {
    if (drctv["containerOptions"] instanceof List) {
      drctv["containerOptions"] = drctv["containerOptions"].join(" ")
    }
    assert drctv["containerOptions"] instanceof CharSequence
  }

  /* DIRECTIVE cpus
    accepted examples:
    - 1
    - 10
  */
  if (drctv.containsKey("cpus")) {
    assert drctv["cpus"] instanceof Integer
  }

  /* DIRECTIVE disk
    accepted examples:
    - "1 GB"
    - "2TB"
    - "3.2KB"
    - "10.B"
  */
  if (drctv.containsKey("disk")) {
    assert drctv["disk"] instanceof CharSequence
    // assert drctv["disk"].matches("[0-9]+(\\.[0-9]*)? *[KMGTPEZY]?B")
    // ^ does not allow closures
  }

  /* DIRECTIVE echo
    accepted examples:
    - true
    - false
  */
  if (drctv.containsKey("echo")) {
    assert drctv["echo"] instanceof Boolean
  }

  /* DIRECTIVE errorStrategy
    accepted examples:
    - "terminate"
    - "finish"
  */
  if (drctv.containsKey("errorStrategy")) {
    assert drctv["errorStrategy"] instanceof CharSequence
    assert drctv["errorStrategy"] in ["terminate", "finish", "ignore", "retry"] : "Unexpected value for errorStrategy"
  }

  /* DIRECTIVE executor
    accepted examples:
    - "local"
    - "sge"
  */
  if (drctv.containsKey("executor")) {
    assert drctv["executor"] instanceof CharSequence
    assert drctv["executor"] in ["local", "sge", "uge", "lsf", "slurm", "pbs", "pbspro", "moab", "condor", "nqsii", "ignite", "k8s", "awsbatch", "google-pipelines"] : "Unexpected value for executor"
  }

  /* DIRECTIVE machineType
    accepted examples:
    - "n1-highmem-8"
  */
  if (drctv.containsKey("machineType")) {
    assert drctv["machineType"] instanceof CharSequence
  }

  /* DIRECTIVE maxErrors
    accepted examples:
    - 1
    - 3
  */
  if (drctv.containsKey("maxErrors")) {
    assert drctv["maxErrors"] instanceof Integer
  }

  /* DIRECTIVE maxForks
    accepted examples:
    - 1
    - 3
  */
  if (drctv.containsKey("maxForks")) {
    assert drctv["maxForks"] instanceof Integer
  }

  /* DIRECTIVE maxRetries
    accepted examples:
    - 1
    - 3
  */
  if (drctv.containsKey("maxRetries")) {
    assert drctv["maxRetries"] instanceof Integer
  }

  /* DIRECTIVE memory
    accepted examples:
    - "1 GB"
    - "2TB"
    - "3.2KB"
    - "10.B"
  */
  if (drctv.containsKey("memory")) {
    assert drctv["memory"] instanceof CharSequence
    // assert drctv["memory"].matches("[0-9]+(\\.[0-9]*)? *[KMGTPEZY]?B")
    // ^ does not allow closures
  }

  /* DIRECTIVE module
    accepted examples:
    - "ncbi-blast/2.2.27"
    - "ncbi-blast/2.2.27:t_coffee/10.0"
    - ["ncbi-blast/2.2.27", "t_coffee/10.0"]
  */
  if (drctv.containsKey("module")) {
    if (drctv["module"] instanceof List) {
      drctv["module"] = drctv["module"].join(":")
    }
    assert drctv["module"] instanceof CharSequence
  }

  /* DIRECTIVE penv
    accepted examples:
    - "smp"
  */
  if (drctv.containsKey("penv")) {
    assert drctv["penv"] instanceof CharSequence
  }

  /* DIRECTIVE pod
    accepted examples:
    - [ label: "key", value: "val" ]
    - [ annotation: "key", value: "val" ]
    - [ env: "key", value: "val" ]
    - [ [label: "l", value: "v"], [env: "e", value: "v"]]
  */
  if (drctv.containsKey("pod")) {
    if (drctv["pod"] instanceof Map) {
      drctv["pod"] = [ drctv["pod"] ]
    }
    assert drctv["pod"] instanceof List
    drctv["pod"].forEach { pod ->
      assert pod instanceof Map
      // TODO: should more checks be added?
      // See https://www.nextflow.io/docs/latest/process.html?highlight=directives#pod
      // e.g. does it contain 'label' and 'value', or 'annotation' and 'value', or ...?
    }
  }

  /* DIRECTIVE publishDir
    accepted examples:
    - []
    - [ [ path: "foo", enabled: true ], [ path: "bar", enabled: false ] ]
    - "/path/to/dir" 
      is transformed to [[ path: "/path/to/dir" ]]
    - [ path: "/path/to/dir", mode: "cache" ]
      is transformed to [[ path: "/path/to/dir", mode: "cache" ]]
  */
  // TODO: should we also look at params["publishDir"]?
  if (drctv.containsKey("publishDir")) {
    def pblsh = drctv["publishDir"]
    
    // check different options
    assert pblsh instanceof List || pblsh instanceof Map || pblsh instanceof CharSequence
    
    // turn into list if not already so
    // for some reason, 'if (!pblsh instanceof List) pblsh = [ pblsh ]' doesn't work.
    pblsh = pblsh instanceof List ? pblsh : [ pblsh ]

    // check elements of publishDir
    pblsh = pblsh.collect{ elem ->
      // turn into map if not already so
      elem = elem instanceof CharSequence ? [ path: elem ] : elem

      // check types and keys
      assert elem instanceof Map : "Expected publish argument '$elem' to be a String or a Map. Found: class ${elem.getClass()}"
      assertMapKeys(elem, [ "path", "mode", "overwrite", "pattern", "saveAs", "enabled" ], ["path"], "publishDir")

      // check elements in map
      assert elem.containsKey("path")
      assert elem["path"] instanceof CharSequence
      if (elem.containsKey("mode")) {
        assert elem["mode"] instanceof CharSequence
        assert elem["mode"] in [ "symlink", "rellink", "link", "copy", "copyNoFollow", "move" ]
      }
      if (elem.containsKey("overwrite")) {
        assert elem["overwrite"] instanceof Boolean
      }
      if (elem.containsKey("pattern")) {
        assert elem["pattern"] instanceof CharSequence
      }
      if (elem.containsKey("saveAs")) {
        assert elem["saveAs"] instanceof CharSequence //: "saveAs as a Closure is currently not supported. Surround your closure with single quotes to get the desired effect. Example: '\{ foo \}'"
      }
      if (elem.containsKey("enabled")) {
        assert elem["enabled"] instanceof Boolean
      }

      // return final result
      elem
    }
    // store final directive
    drctv["publishDir"] = pblsh
  }

  /* DIRECTIVE queue
    accepted examples:
    - "long"
    - "short,long"
    - ["short", "long"]
  */
  if (drctv.containsKey("queue")) {
    if (drctv["queue"] instanceof List) {
      drctv["queue"] = drctv["queue"].join(",")
    }
    assert drctv["queue"] instanceof CharSequence
  }

  /* DIRECTIVE label
    accepted examples:
    - "big_mem"
    - "big_cpu"
    - ["big_mem", "big_cpu"]
  */
  if (drctv.containsKey("label")) {
    if (drctv["label"] instanceof CharSequence) {
      drctv["label"] = [ drctv["label"] ]
    }
    assert drctv["label"] instanceof List
    drctv["label"].forEach { label ->
      assert label instanceof CharSequence
      // assert label.matches("[a-zA-Z0-9]([a-zA-Z0-9_]*[a-zA-Z0-9])?")
      // ^ does not allow closures
    }
  }

  /* DIRECTIVE scratch
    accepted examples:
    - true
    - "/path/to/scratch"
    - '$MY_PATH_TO_SCRATCH'
    - "ram-disk"
  */
  if (drctv.containsKey("scratch")) {
    assert drctv["scratch"] == true || drctv["scratch"] instanceof CharSequence
  }

  /* DIRECTIVE storeDir
    accepted examples:
    - "/path/to/storeDir"
  */
  if (drctv.containsKey("storeDir")) {
    assert drctv["storeDir"] instanceof CharSequence
  }

  /* DIRECTIVE stageInMode
    accepted examples:
    - "copy"
    - "link"
  */
  if (drctv.containsKey("stageInMode")) {
    assert drctv["stageInMode"] instanceof CharSequence
    assert drctv["stageInMode"] in ["copy", "link", "symlink", "rellink"]
  }

  /* DIRECTIVE stageOutMode
    accepted examples:
    - "copy"
    - "link"
  */
  if (drctv.containsKey("stageOutMode")) {
    assert drctv["stageOutMode"] instanceof CharSequence
    assert drctv["stageOutMode"] in ["copy", "move", "rsync"]
  }

  /* DIRECTIVE tag
    accepted examples:
    - "foo"
    - '$id'
  */
  if (drctv.containsKey("tag")) {
    assert drctv["tag"] instanceof CharSequence
  }

  /* DIRECTIVE time
    accepted examples:
    - "1h"
    - "2days"
    - "1day 6hours 3minutes 30seconds"
  */
  if (drctv.containsKey("time")) {
    assert drctv["time"] instanceof CharSequence
    // todo: validation regex?
  }

  return drctv
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/processWorkflowArgs.nf'
// depends on: thisConfig, thisDefaultWorkflowArgs
def processWorkflowArgs(Map args) {
  // override defaults with args
  def workflowArgs = thisDefaultWorkflowArgs + args

  // check whether 'key' exists
  assert workflowArgs.containsKey("key") : "Error in module '${thisConfig.functionality.name}': key is a required argument"

  // if 'key' is a closure, apply it to the original key
  if (workflowArgs["key"] instanceof Closure) {
    workflowArgs["key"] = workflowArgs["key"](thisConfig.functionality.name)
  }
  def key = workflowArgs["key"]
  assert key instanceof CharSequence : "Expected process argument 'key' to be a String. Found: class ${key.getClass()}"
  assert key ==~ /^[a-zA-Z_]\w*$/ : "Error in module '$key': Expected process argument 'key' to consist of only letters, digits or underscores. Found: ${key}"

  // check for any unexpected keys
  def expectedKeys = ["key", "directives", "auto", "map", "mapId", "mapData", "mapPassthrough", "filter", "fromState", "toState", "args", "renameKeys", "debug"]
  def unexpectedKeys = workflowArgs.keySet() - expectedKeys
  assert unexpectedKeys.isEmpty() : "Error in module '$key': unexpected arguments to the '.run()' function: '${unexpectedKeys.join("', '")}'"

  // check whether directives exists and apply defaults
  assert workflowArgs.containsKey("directives") : "Error in module '$key': directives is a required argument"
  assert workflowArgs["directives"] instanceof Map : "Error in module '$key': Expected process argument 'directives' to be a Map. Found: class ${workflowArgs['directives'].getClass()}"
  workflowArgs["directives"] = processDirectives(thisDefaultWorkflowArgs.directives + workflowArgs["directives"])

  // check whether directives exists and apply defaults
  assert workflowArgs.containsKey("auto") : "Error in module '$key': auto is a required argument"
  assert workflowArgs["auto"] instanceof Map : "Error in module '$key': Expected process argument 'auto' to be a Map. Found: class ${workflowArgs['auto'].getClass()}"
  workflowArgs["auto"] = processAuto(thisDefaultWorkflowArgs.auto + workflowArgs["auto"])

  // auto define publish, if so desired
  if (workflowArgs.auto.publish == true && (workflowArgs.directives.publishDir != null ? workflowArgs.directives.publishDir : [:]).isEmpty()) {
    // can't assert at this level thanks to the no_publish profile
    // assert params.containsKey("publishDir") || params.containsKey("publish_dir") : 
    //   "Error in module '${workflowArgs['key']}': if auto.publish is true, params.publish_dir needs to be defined.\n" +
    //   "  Example: params.publish_dir = \"./output/\""
    def publishDir = getPublishDir()
    
    if (publishDir != null) {
      workflowArgs.directives.publishDir = [[ 
        path: publishDir, 
        saveAs: "{ it.startsWith('.') ? null : it }", // don't publish hidden files, by default
        mode: "copy"
      ]]
    }
  }

  // auto define transcript, if so desired
  if (workflowArgs.auto.transcript == true) {
    // can't assert at this level thanks to the no_publish profile
    // assert params.containsKey("transcriptsDir") || params.containsKey("transcripts_dir") || params.containsKey("publishDir") || params.containsKey("publish_dir") : 
    //   "Error in module '${workflowArgs['key']}': if auto.transcript is true, either params.transcripts_dir or params.publish_dir needs to be defined.\n" +
    //   "  Example: params.transcripts_dir = \"./transcripts/\""
    def transcriptsDir = 
      params.containsKey("transcripts_dir") ? params.transcripts_dir : 
      params.containsKey("transcriptsDir") ? params.transcriptsDir : 
      params.containsKey("publish_dir") ? params.publish_dir + "/_transcripts" :
      params.containsKey("publishDir") ? params.publishDir + "/_transcripts" : 
      null
    if (transcriptsDir != null) {
      def timestamp = nextflow.Nextflow.getSession().getWorkflowMetadata().start.format('yyyy-MM-dd_HH-mm-ss')
      def transcriptsPublishDir = [ 
        path: "$transcriptsDir/$timestamp/\${task.process.replaceAll(':', '-')}/\${id}/",
        saveAs: "{ it.startsWith('.') ? it.replaceAll('^.', '') : null }", 
        mode: "copy"
      ]
      def publishDirs = workflowArgs.directives.publishDir != null ? workflowArgs.directives.publishDir : null ? workflowArgs.directives.publishDir : []
      workflowArgs.directives.publishDir = publishDirs + transcriptsPublishDir
    }
  }

  // if this is a stubrun, remove certain directives?
  if (workflow.stubRun) {
    workflowArgs.directives.keySet().removeAll(["publishDir", "cpus", "memory", "label"])
  }

  for (nam in ["map", "mapId", "mapData", "mapPassthrough", "filter"]) {
    if (workflowArgs.containsKey(nam) && workflowArgs[nam]) {
      assert workflowArgs[nam] instanceof Closure : "Error in module '$key': Expected process argument '$nam' to be null or a Closure. Found: class ${workflowArgs[nam].getClass()}"
    }
  }

  // TODO: should functions like 'map', 'mapId', 'mapData', 'mapPassthrough' be deprecated as well?
  for (nam in ["renameKeys"]) {
    if (workflowArgs.containsKey(nam) && workflowArgs[nam] != null) {
      log.warn "module '$key': workflow argument '$nam' will be deprecated in Viash 0.9.0. Please use 'fromState' and 'toState' instead."
    }
  }

  // check fromState
  workflowArgs["fromState"] = _processFromState(workflowArgs.get("fromState"), key, thisConfig)

  // check toState
  workflowArgs["toState"] = _processToState(workflowArgs.get("toState"), key, thisConfig)

  // return output
  return workflowArgs
}

def _processFromState(fromState, key_, config_) {
  assert fromState == null || fromState instanceof Closure || fromState instanceof Map || fromState instanceof List :
    "Error in module '$key_': Expected process argument 'fromState' to be null, a Closure, a Map, or a List. Found: class ${fromState.getClass()}"
  if (fromState == null) {
    return null
  }
  
  // if fromState is a List, convert to map
  if (fromState instanceof List) {
    // check whether fromstate is a list[string]
    assert fromState.every{it instanceof CharSequence} : "Error in module '$key_': fromState is a List, but not all elements are Strings"
    fromState = fromState.collectEntries{[it, it]}
  }

  // if fromState is a map, convert to closure
  if (fromState instanceof Map) {
    // check whether fromstate is a map[string, string]
    assert fromState.values().every{it instanceof CharSequence} : "Error in module '$key_': fromState is a Map, but not all values are Strings"
    assert fromState.keySet().every{it instanceof CharSequence} : "Error in module '$key_': fromState is a Map, but not all keys are Strings"
    def fromStateMap = fromState.clone()
    def requiredInputNames = thisConfig.functionality.allArguments.findAll{it.required && it.direction == "Input"}.collect{it.plainName}
    // turn the map into a closure to be used later on
    fromState = { it ->
      def state = it[1]
      assert state instanceof Map : "Error in module '$key_': the state is not a Map"
      def data = fromStateMap.collectMany{newkey, origkey ->
        // check whether newkey corresponds to a required argument
        if (state.containsKey(origkey)) {
          [[newkey, state[origkey]]]
        } else if (!requiredInputNames.contains(origkey)) {
          []
        } else {
          throw new Exception("Error in module '$key_': fromState key '$origkey' not found in current state")
        }
      }.collectEntries()
      data
    }
  }
  
  return fromState
}

def _processToState(toState, key_, config_) {
  if (toState == null) {
    toState = { tup -> tup[1] }
  }

  // toState should be a closure, map[string, string], or list[string]
  assert toState instanceof Closure || toState instanceof Map || toState instanceof List :
    "Error in module '$key_': Expected process argument 'toState' to be a Closure, a Map, or a List. Found: class ${toState.getClass()}"

  // if toState is a List, convert to map
  if (toState instanceof List) {
    // check whether toState is a list[string]
    assert toState.every{it instanceof CharSequence} : "Error in module '$key_': toState is a List, but not all elements are Strings"
    toState = toState.collectEntries{[it, it]}
  }

  // if toState is a map, convert to closure
  if (toState instanceof Map) {
    // check whether toState is a map[string, string]
    assert toState.values().every{it instanceof CharSequence} : "Error in module '$key_': toState is a Map, but not all values are Strings"
    assert toState.keySet().every{it instanceof CharSequence} : "Error in module '$key_': toState is a Map, but not all keys are Strings"
    def toStateMap = toState.clone()
    def requiredOutputNames = config_.functionality.allArguments.findAll{it.required && it.direction == "Output"}.collect{it.plainName}
    // turn the map into a closure to be used later on
    toState = { it ->
      def output = it[1]
      def state = it[2]
      assert output instanceof Map : "Error in module '$key_': the output is not a Map"
      assert state instanceof Map : "Error in module '$key_': the state is not a Map"
      def extraEntries = toStateMap.collectMany{newkey, origkey ->
        // check whether newkey corresponds to a required argument
        if (output.containsKey(origkey)) {
          [[newkey, output[origkey]]]
        } else if (!requiredOutputNames.contains(origkey)) {
          []
        } else {
          throw new Exception("Error in module '$key_': toState key '$origkey' not found in current output")
        }
      }.collectEntries()
      state + extraEntries
    }
  }

  return toState
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/vdsl3ProcessFactory.nf'
// depends on: thisConfig, thisScript, session?
def vdsl3ProcessFactory(Map workflowArgs) {
  // autodetect process key
  def wfKey = workflowArgs["key"]
  def procKeyPrefix = "${wfKey}_process"
  def meta = nextflow.script.ScriptMeta.current()
  def existing = meta.getProcessNames().findAll{it.startsWith(procKeyPrefix)}
  def numbers = existing.collect{it.replace(procKeyPrefix, "0").toInteger()}
  def newNumber = (numbers + [-1]).max() + 1

  def procKey = newNumber == 0 ? procKeyPrefix : "$procKeyPrefix$newNumber"

  if (newNumber > 0) {
    log.warn "Key for module '${wfKey}' is duplicated.\n",
      "If you run a component multiple times in the same workflow,\n" +
      "it's recommended you set a unique key for every call,\n" +
      "for example: ${wfKey}.run(key: \"foo\")."
  }

  // subset directives and convert to list of tuples
  def drctv = workflowArgs.directives

  // TODO: unit test the two commands below
  // convert publish array into tags
  def valueToStr = { val ->
    // ignore closures
    if (val instanceof CharSequence) {
      if (!val.matches('^[{].*[}]$')) {
        '"' + val + '"'
      } else {
        val
      }
    } else if (val instanceof List) {
      "[" + val.collect{valueToStr(it)}.join(", ") + "]"
    } else if (val instanceof Map) {
      "[" + val.collect{k, v -> k + ": " + valueToStr(v)}.join(", ") + "]"
    } else {
      val.inspect()
    }
  }

  // multiple entries allowed: label, publishdir
  def drctvStrs = drctv.collect { key, value ->
    if (key in ["label", "publishDir"]) {
      value.collect{ val ->
        if (val instanceof Map) {
          "\n$key " + val.collect{ k, v -> k + ": " + valueToStr(v) }.join(", ")
        } else if (val == null) {
          ""
        } else {
          "\n$key " + valueToStr(val)
        }
      }.join()
    } else if (value instanceof Map) {
      "\n$key " + value.collect{ k, v -> k + ": " + valueToStr(v) }.join(", ")
    } else {
      "\n$key " + valueToStr(value)
    }
  }.join()

  def inputPaths = thisConfig.functionality.allArguments
    .findAll { it.type == "file" && it.direction == "input" }
    .collect { ', path(viash_par_' + it.plainName + ', stageAs: "_viash_par/' + it.plainName + '_?/*")' }
    .join()

  def outputPaths = thisConfig.functionality.allArguments
    .findAll { it.type == "file" && it.direction == "output" }
    .collect { par ->
      // insert dummy into every output (see nextflow-io/nextflow#2678)
      if (!par.multiple) {
        ', path{[".exitcode", args.' + par.plainName + ']}'
      } else {
        ', path{[".exitcode"] + args.' + par.plainName + '}'
      }
    }
    .join()

  // TODO: move this functionality somewhere else?
  if (workflowArgs.auto.transcript) {
    outputPaths = outputPaths + ', path{[".exitcode", ".command*"]}'
  } else {
    outputPaths = outputPaths + ', path{[".exitcode"]}'
  }

  // create dirs for output files (based on BashWrapper.createParentFiles)
  def createParentStr = thisConfig.functionality.allArguments
    .findAll { it.type == "file" && it.direction == "output" && it.create_parent }
    .collect { par -> 
      "\${ args.containsKey(\"${par.plainName}\") ? \"mkdir_parent \\\"\" + (args[\"${par.plainName}\"] instanceof String ? args[\"${par.plainName}\"] : args[\"${par.plainName}\"].join('\" \"')) + \"\\\"\" : \"\" }"
    }
    .join("\n")

  // construct inputFileExports
  def inputFileExports = thisConfig.functionality.allArguments
    .findAll { it.type == "file" && it.direction.toLowerCase() == "input" }
    .collect { par ->
      viash_par_contents = "(viash_par_${par.plainName} instanceof List ? viash_par_${par.plainName}.join(\"${par.multiple_sep}\") : viash_par_${par.plainName})"
      "\n\${viash_par_${par.plainName}.empty ? \"\" : \"export VIASH_PAR_${par.plainName.toUpperCase()}=\\\"\" + ${viash_par_contents} + \"\\\"\"}"
    }

  // NOTE: if using docker, use /tmp instead of tmpDir!
  def tmpDir = java.nio.file.Paths.get(
    System.getenv('NXF_TEMP') ?: 
    System.getenv('VIASH_TEMP') ?: 
    System.getenv('VIASH_TMPDIR') ?: 
    System.getenv('VIASH_TEMPDIR') ?: 
    System.getenv('VIASH_TMP') ?: 
    System.getenv('TEMP') ?: 
    System.getenv('TMPDIR') ?: 
    System.getenv('TEMPDIR') ?:
    System.getenv('TMP') ?: 
    '/tmp'
  ).toAbsolutePath()

  // construct stub
  def stub = thisConfig.functionality.allArguments
    .findAll { it.type == "file" && it.direction == "output" }
    .collect { par -> 
      "\${ args.containsKey(\"${par.plainName}\") ? \"touch2 \\\"\" + (args[\"${par.plainName}\"] instanceof String ? args[\"${par.plainName}\"].replace(\"_*\", \"_0\") : args[\"${par.plainName}\"].join('\" \"')) + \"\\\"\" : \"\" }"
    }
    .join("\n")

  // escape script
  def escapedScript = thisScript.replace('\\', '\\\\').replace('$', '\\$').replace('"""', '\\"\\"\\"')

  // publishdir assert
  def assertStr = (workflowArgs.auto.publish == true) || workflowArgs.auto.transcript ? 
    """\nassert task.publishDir.size() > 0: "if auto.publish is true, params.publish_dir needs to be defined.\\n  Example: --publish_dir './output/'" """ :
    ""

  // generate process string
  def procStr = 
  """nextflow.enable.dsl=2
  |
  |process $procKey {$drctvStrs
  |input:
  |  tuple val(id)$inputPaths, val(args), path(resourcesDir, stageAs: ".viash_meta_resources")
  |output:
  |  tuple val("\$id")$outputPaths, optional: true
  |stub:
  |\"\"\"
  |touch2() { mkdir -p "\\\$(dirname "\\\$1")" && touch "\\\$1" ; }
  |$stub
  |\"\"\"
  |script:$assertStr
  |def escapeText = { s -> s.toString().replaceAll('([`"])', '\\\\\\\\\$1') }
  |def parInject = args
  |  .findAll{key, value -> value != null}
  |  .collect{key, value -> "export VIASH_PAR_\${key.toUpperCase()}=\\\"\${escapeText(value)}\\\""}
  |  .join("\\n")
  |\"\"\"
  |# meta exports
  |# export VIASH_META_RESOURCES_DIR="\${resourcesDir.toRealPath().toAbsolutePath()}"
  |export VIASH_META_RESOURCES_DIR="\${resourcesDir}"
  |export VIASH_META_TEMP_DIR="${['docker', 'podman', 'charliecloud'].any{ it == workflow.containerEngine } ? '/tmp' : tmpDir}"
  |export VIASH_META_FUNCTIONALITY_NAME="${thisConfig.functionality.name}"
  |# export VIASH_META_EXECUTABLE="\\\$VIASH_META_RESOURCES_DIR/\\\$VIASH_META_FUNCTIONALITY_NAME"
  |export VIASH_META_CONFIG="\\\$VIASH_META_RESOURCES_DIR/.config.vsh.yaml"
  |\${task.cpus ? "export VIASH_META_CPUS=\$task.cpus" : "" }
  |\${task.memory?.bytes != null ? "export VIASH_META_MEMORY_B=\$task.memory.bytes" : "" }
  |if [ ! -z \\\${VIASH_META_MEMORY_B+x} ]; then
  |  export VIASH_META_MEMORY_KB=\\\$(( (\\\$VIASH_META_MEMORY_B+1023) / 1024 ))
  |  export VIASH_META_MEMORY_MB=\\\$(( (\\\$VIASH_META_MEMORY_KB+1023) / 1024 ))
  |  export VIASH_META_MEMORY_GB=\\\$(( (\\\$VIASH_META_MEMORY_MB+1023) / 1024 ))
  |  export VIASH_META_MEMORY_TB=\\\$(( (\\\$VIASH_META_MEMORY_GB+1023) / 1024 ))
  |  export VIASH_META_MEMORY_PB=\\\$(( (\\\$VIASH_META_MEMORY_TB+1023) / 1024 ))
  |fi
  |
  |# meta synonyms
  |export VIASH_TEMP="\\\$VIASH_META_TEMP_DIR"
  |export TEMP_DIR="\\\$VIASH_META_TEMP_DIR"
  |
  |# create output dirs if need be
  |function mkdir_parent {
  |  for file in "\\\$@"; do 
  |    mkdir -p "\\\$(dirname "\\\$file")"
  |  done
  |}
  |$createParentStr
  |
  |# argument exports${inputFileExports.join()}
  |\$parInject
  |
  |# process script
  |${escapedScript}
  |\"\"\"
  |}
  |""".stripMargin()

  // TODO: print on debug
  // if (workflowArgs.debug == true) {
  //   println("######################\n$procStr\n######################")
  // }

  // create runtime process
  def ownerParams = new nextflow.script.ScriptBinding.ParamsMap()
  def binding = new nextflow.script.ScriptBinding().setParams(ownerParams)
  def module = new nextflow.script.IncludeDef.Module(name: procKey)
  def scriptParser = new nextflow.script.ScriptParser(session)
    .setModule(true)
    .setBinding(binding)
  scriptParser.scriptPath = meta.getScriptPath()
  def moduleScript = scriptParser.runScript(procStr)
    .getScript()

  // register module in meta
  meta.addModule(moduleScript, module.name, module.alias)

  // retrieve and return process from meta
  return meta.getProcess(procKey)
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/vdsl3WorkflowFactory.nf'
// depends on: thisConfig, resourcesDir
def vdsl3WorkflowFactory(Map args) {
  def key = args["key"]
  def processObj = null

  workflow processWf {
    take: input_
    main:

    if (processObj == null) {
      processObj = vdsl3ProcessFactory(args)
    }
    
    output_ = input_
      | map { tuple ->
        def id = tuple[0]
        def data_ = tuple[1]

        if (workflow.stubRun) {
          // add id if missing
          data_ = [id: 'stub'] + data_
        }

        // process input files separately
        def inputPaths = thisConfig.functionality.allArguments
          .findAll { it.type == "file" && it.direction == "input" }
          .collect { par ->
            def val = data_.containsKey(par.plainName) ? data_[par.plainName] : []
            def inputFiles = []
            if (val == null) {
              inputFiles = []
            } else if (val instanceof List) {
              inputFiles = val
            } else if (val instanceof Path) {
              inputFiles = [ val ]
            } else {
              inputFiles = []
            }
            if (!workflow.stubRun) {
              // throw error when an input file doesn't exist
              inputFiles.each{ file -> 
                assert file.exists() :
                  "Error in module '${key}' id '${id}' argument '${par.plainName}'.\n" +
                  "  Required input file does not exist.\n" +
                  "  Path: '$file'.\n" +
                  "  Expected input file to exist"
              }
            }
            inputFiles 
          } 

        // remove input files
        def argsExclInputFiles = thisConfig.functionality.allArguments
          .findAll { (it.type != "file" || it.direction != "input") && data_.containsKey(it.plainName) }
          .collectEntries { par ->
            def parName = par.plainName
            def val = data_[parName]
            if (par.multiple && val instanceof Collection) {
              val = val.join(par.multiple_sep)
            }
            if (par.direction == "output" && par.type == "file") {
              val = val.replaceAll('\\$id', id).replaceAll('\\$key', key)
            }
            [parName, val]
          }

        [ id ] + inputPaths + [ argsExclInputFiles, resourcesDir ]
      }
      | processObj
      | map { output ->
        def outputFiles = thisConfig.functionality.allArguments
          .findAll { it.type == "file" && it.direction == "output" }
          .indexed()
          .collectEntries{ index, par ->
            out = output[index + 1]
            // strip dummy '.exitcode' file from output (see nextflow-io/nextflow#2678)
            if (!out instanceof List || out.size() <= 1) {
              if (par.multiple) {
                out = []
              } else {
                assert !par.required :
                    "Error in module '${key}' id '${output[0]}' argument '${par.plainName}'.\n" +
                    "  Required output file is missing"
                out = null
              }
            } else if (out.size() == 2 && !par.multiple) {
              out = out[1]
            } else {
              out = out.drop(1)
            }
            [ par.plainName, out ]
          }
        
        // drop null outputs
        outputFiles.removeAll{it.value == null}

        [ output[0], outputFiles ]
      }
    emit: output_
  }

  return processWf
}

// helper file: 'src/main/resources/io/viash/platforms/nextflow/workflowFactory/workflowFactory.nf'
def _debug(workflowArgs, debugKey) {
  if (workflowArgs.debug) {
    view { "process '${workflowArgs.key}' $debugKey tuple: $it"  }
  } else {
    map { it }
  }
}

// depends on: thisConfig, innerWorkflowFactory
def workflowFactory(Map args) {
  def workflowArgs = processWorkflowArgs(args)
  def key_ = workflowArgs["key"]
  
  workflow workflowInstance {
    take: input_

    main:
    mid1_ = input_
      | checkUniqueIds([:])
      | _debug(workflowArgs, "input")
      | map { tuple ->
        tuple = tuple.clone()
        
        if (workflowArgs.map) {
          tuple = workflowArgs.map(tuple)
        }
        if (workflowArgs.mapId) {
          tuple[0] = workflowArgs.mapId(tuple[0])
        }
        if (workflowArgs.mapData) {
          tuple[1] = workflowArgs.mapData(tuple[1])
        }
        if (workflowArgs.mapPassthrough) {
          tuple = tuple.take(2) + workflowArgs.mapPassthrough(tuple.drop(2))
        }

        // check tuple
        assert tuple instanceof List : 
          "Error in module '${key_}': element in channel should be a tuple [id, data, ...otherargs...]\n" +
          "  Example: [\"id\", [input: file('foo.txt'), arg: 10]].\n" +
          "  Expected class: List. Found: tuple.getClass() is ${tuple.getClass()}"
        assert tuple.size() >= 2 : 
          "Error in module '${key_}': expected length of tuple in input channel to be two or greater.\n" +
          "  Example: [\"id\", [input: file('foo.txt'), arg: 10]].\n" +
          "  Found: tuple.size() == ${tuple.size()}"
        
        // check id field
        assert tuple[0] instanceof CharSequence : 
          "Error in module '${key_}': first element of tuple in channel should be a String\n" +
          "  Example: [\"id\", [input: file('foo.txt'), arg: 10]].\n" +
          "  Found: ${tuple[0]}"
        
        // match file to input file
        if (workflowArgs.auto.simplifyInput && (tuple[1] instanceof Path || tuple[1] instanceof List)) {
          def inputFiles = thisConfig.functionality.allArguments
            .findAll { it.type == "file" && it.direction == "input" }
          
          assert inputFiles.size() == 1 : 
              "Error in module '${key_}' id '${tuple[0]}'.\n" +
              "  Anonymous file inputs are only allowed when the process has exactly one file input.\n" +
              "  Expected: inputFiles.size() == 1. Found: inputFiles.size() is ${inputFiles.size()}"

          tuple[1] = [[ inputFiles[0].plainName, tuple[1] ]].collectEntries()
        }

        // check data field
        assert tuple[1] instanceof Map : 
          "Error in module '${key_}' id '${tuple[0]}': second element of tuple in channel should be a Map\n" +
          "  Example: [\"id\", [input: file('foo.txt'), arg: 10]].\n" +
          "  Expected class: Map. Found: tuple[1].getClass() is ${tuple[1].getClass()}"

        // rename keys of data field in tuple
        if (workflowArgs.renameKeys) {
          assert workflowArgs.renameKeys instanceof Map : 
              "Error renaming data keys in module '${key_}' id '${tuple[0]}'.\n" +
              "  Example: renameKeys: ['new_key': 'old_key'].\n" +
              "  Expected class: Map. Found: renameKeys.getClass() is ${workflowArgs.renameKeys.getClass()}"
          assert tuple[1] instanceof Map : 
              "Error renaming data keys in module '${key_}' id '${tuple[0]}'.\n" +
              "  Expected class: Map. Found: tuple[1].getClass() is ${tuple[1].getClass()}"

          // TODO: allow renameKeys to be a function?
          workflowArgs.renameKeys.each { newKey, oldKey ->
            assert newKey instanceof CharSequence : 
              "Error renaming data keys in module '${key_}' id '${tuple[0]}'.\n" +
              "  Example: renameKeys: ['new_key': 'old_key'].\n" +
              "  Expected class of newKey: String. Found: newKey.getClass() is ${newKey.getClass()}"
            assert oldKey instanceof CharSequence : 
              "Error renaming data keys in module '${key_}' id '${tuple[0]}'.\n" +
              "  Example: renameKeys: ['new_key': 'old_key'].\n" +
              "  Expected class of oldKey: String. Found: oldKey.getClass() is ${oldKey.getClass()}"
            assert tuple[1].containsKey(oldKey) : 
              "Error renaming data keys in module '${key}' id '${tuple[0]}'.\n" +
              "  Key '$oldKey' is missing in the data map. tuple[1].keySet() is '${tuple[1].keySet()}'"
            tuple[1].put(newKey, tuple[1][oldKey])
          }
          tuple[1].keySet().removeAll(workflowArgs.renameKeys.collect{ newKey, oldKey -> oldKey })
        }
        tuple
      }

    if (workflowArgs.filter) {
      mid2_ = mid1_
        | filter{workflowArgs.filter(it)}
    } else {
      mid2_ = mid1_
    }

    if (workflowArgs.fromState) {
      mid3_ = mid2_
        | map{
          def new_data = workflowArgs.fromState(it.take(2))
          [it[0], new_data]
        }
    } else {
      mid3_ = mid2_
    }

    // fill in defaults
    mid4_ = mid3_
      | map { tuple ->
        def id_ = tuple[0]
        def data_ = tuple[1]

        // TODO: could move fromState to here

        // fetch default params from functionality
        def defaultArgs = thisConfig.functionality.allArguments
          .findAll { it.containsKey("default") }
          .collectEntries { [ it.plainName, it.default ] }

        // fetch overrides in params
        def paramArgs = thisConfig.functionality.allArguments
          .findAll { par ->
            def argKey = key_ + "__" + par.plainName
            params.containsKey(argKey)
          }
          .collectEntries { [ it.plainName, params[key_ + "__" + it.plainName] ] }
        
        // fetch overrides in data
        def dataArgs = thisConfig.functionality.allArguments
          .findAll { data_.containsKey(it.plainName) }
          .collectEntries { [ it.plainName, data_[it.plainName] ] }
        
        // combine params
        def combinedArgs = defaultArgs + paramArgs + workflowArgs.args + dataArgs

        // remove arguments with explicit null values
        combinedArgs
          .removeAll{_, val -> val == null || val == "viash_no_value" || val == "force_null"}

        combinedArgs = processInputs(combinedArgs, thisConfig, id_, key_)

        [id_, combinedArgs] + tuple.drop(2)
      }

    // TODO: move some of the _meta.join_id wrangling to the safeJoin() function.

    out0_ = mid4_
      | _debug(workflowArgs, "processed")
      // run workflow
      | innerWorkflowFactory(workflowArgs)
      // check output tuple
      | map { id_, output_ ->

        // see if output map contains metadata
        def meta_ =
          output_ instanceof Map && output_.containsKey("_meta") ? 
          output_["_meta"] :
          [:]
        if (!meta_.containsKey("join_id")) {
          meta_ = meta_ + ["join_id": id_]
        }
        
        // remove metadata
        output_ = output_.findAll{k, v -> k != "_meta"}

        output_ = processOutputs(output_, thisConfig, id_, key_)

        if (workflowArgs.auto.simplifyOutput && output_.size() == 1) {
          output_ = output_.values()[0]
        }

        [meta_.join_id, meta_, id_, output_]
      }
      // | view{"out0_: ${it.take(3)}"}

    // TODO: this join will fail if the keys changed during the innerWorkflowFactory
    // join the output [join_id, meta, id, output] with the previous state [id, state, ...]
    out1_ = safeJoin(out0_, mid2_, key_)
      // input tuple format: [join_id, meta, id, output, prev_state, ...]
      // output tuple format: [join_id, meta, id, new_state, ...]
      | map{ tup ->
        def new_state = workflowArgs.toState(tup.drop(2).take(3))
        tup.take(3) + [new_state] + tup.drop(5)
      }

    if (workflowArgs.auto.publish == "state") {
      out1pub_ = out1_
        // input tuple format: [join_id, meta, id, new_state, ...]
        // output tuple format: [join_id, meta, id, new_state]
        | map{ tup ->
          tup.take(4)
        }

      safeJoin(out1pub_, mid4_, key_)
        // input tuple format: [join_id, meta, id, new_state, orig_state, ...]
        // output tuple format: [id, new_state, orig_state]
        | map { tup ->
          tup.drop(2).take(3)
      }
        | publishStatesByConfig(key: key_, config: thisConfig)
    }

    // remove join_id and meta
    out2_ = out1_
      | map { tup ->
        // input tuple format: [join_id, meta, id, new_state, ...]
        // output tuple format: [id, new_state, ...]
        tup.drop(2)
      }
      | _debug(workflowArgs, "output")

    out2_

    emit: out2_
  }

  def wf = workflowInstance.cloneWithName(key_)

  // add factory function
  wf.metaClass.run = { runArgs ->
    workflowFactory(runArgs)
  }
  // add config to module for later introspection
  wf.metaClass.config = thisConfig

  return wf
}
