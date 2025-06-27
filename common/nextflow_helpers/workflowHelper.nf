////////////////////////////
// VDSL3 helper functions //
////////////////////////////

// helper file: 'src/main/resources/io/viash/runners/nextflow/arguments/_checkArgumentType.nf'
class UnexpectedArgumentTypeException extends Exception {
  String errorIdentifier
  String stage
  String plainName
  String expectedClass
  String foundClass
  
  // ${key ? " in module '$key'" : ""}${id ? " id '$id'" : ""}
  UnexpectedArgumentTypeException(String errorIdentifier, String stage, String plainName, String expectedClass, String foundClass) {
    super("Error${errorIdentifier ? " $errorIdentifier" : ""}:${stage ? " $stage" : "" } argument '${plainName}' has the wrong type. " +
      "Expected type: ${expectedClass}. Found type: ${foundClass}")
    this.errorIdentifier = errorIdentifier
    this.stage = stage
    this.plainName = plainName
    this.expectedClass = expectedClass
    this.foundClass = foundClass
  }
}

/**
  * Checks if the given value is of the expected type. If not, an exception is thrown.
  *
  * @param stage The stage of the argument (input or output)
  * @param par The parameter definition
  * @param value The value to check
  * @param errorIdentifier The identifier to use in the error message
  * @return The value, if it is of the expected type
  * @throws UnexpectedArgumentTypeException If the value is not of the expected type
*/
def _checkArgumentType(String stage, Map par, Object value, String errorIdentifier) {
  // expectedClass will only be != null if value is not of the expected type
  def expectedClass = null
  def foundClass = null
  
  // todo: split if need be
  
  if (!par.required && value == null) {
    expectedClass = null
  } else if (par.multiple) {
    if (value !instanceof Collection) {
      value = [value]
    }
    
    // split strings
    value = value.collectMany{ val ->
      if (val instanceof String) {
        // collect() to ensure that the result is a List and not simply an array
        val.split(par.multiple_sep).collect()
      } else {
        [val]
      }
    }

    // process globs
    if (par.type == "file" && par.direction == "input") {
      value = value.collect{ it instanceof String ? file(it, hidden: true) : it }.flatten()
    }

    // check types of elements in list
    try {
      value = value.collect { listVal ->
        _checkArgumentType(stage, par + [multiple: false], listVal, errorIdentifier)
      }
    } catch (UnexpectedArgumentTypeException e) {
      expectedClass = "List[${e.expectedClass}]"
      foundClass = "List[${e.foundClass}]"
    }
  } else if (par.type == "string") {
    // cast to string if need be
    if (value instanceof GString) {
      value = value.toString()
    }
    expectedClass = value instanceof String ? null : "String"
  } else if (par.type == "integer") {
    // cast to integer if need be
    if (value instanceof String) {
      try {
        value = value.toInteger()
      } catch (NumberFormatException e) {
        // do nothing
      }
    }
    if (value instanceof java.math.BigInteger) {
      value = value.intValue()
    }
    expectedClass = value instanceof Integer ? null : "Integer"
  } else if (par.type == "long") {
    // cast to long if need be
    if (value instanceof String) {
      try {
        value = value.toLong()
      } catch (NumberFormatException e) {
        // do nothing
      }
    }
    if (value instanceof Integer) {
      value = value.toLong()
    }
    expectedClass = value instanceof Long ? null : "Long"
  } else if (par.type == "double") {
    // cast to double if need be
    if (value instanceof String) {
      try {
        value = value.toDouble()
      } catch (NumberFormatException e) {
        // do nothing
      }
    }
    if (value instanceof java.math.BigDecimal) {
      value = value.doubleValue()
    }
    if (value instanceof Float) {
      value = value.toDouble()
    }
    expectedClass = value instanceof Double ? null : "Double"
  } else if (par.type == "boolean" | par.type == "boolean_true" | par.type == "boolean_false") {
    // cast to boolean if need be
    if (value instanceof String) {
      def valueLower = value.toLowerCase()
      if (valueLower == "true") {
        value = true
      } else if (valueLower == "false") {
        value = false
      }
    }
    expectedClass = value instanceof Boolean ? null : "Boolean"
  } else if (par.type == "file" && (par.direction == "input" || stage == "output")) {
    // cast to path if need be
    if (value instanceof String) {
      value = file(value, hidden: true)
    }
    if (value instanceof File) {
      value = value.toPath()
    }
    expectedClass = value instanceof Path ? null : "Path"
  } else if (par.type == "file" && stage == "input" && par.direction == "output") {
    // cast to string if need be
    if (value instanceof GString) {
      value = value.toString()
    }
    expectedClass = value instanceof String ? null : "String"
  } else {
    // didn't find a match for par.type
    expectedClass = par.type
  }

  if (expectedClass != null) {
    if (foundClass == null) {
      foundClass = value.getClass().getName()
    }
    throw new UnexpectedArgumentTypeException(errorIdentifier, stage, par.plainName, expectedClass, foundClass)
  }
  
  return value
}
// helper file: 'src/main/resources/io/viash/runners/nextflow/arguments/_processInputValues.nf'
Map _processInputValues(Map inputs, Map config, String id, String key) {
  if (!workflow.stubRun) {
    config.allArguments.each { arg ->
      if (arg.required) {
        assert inputs.containsKey(arg.plainName) && inputs.get(arg.plainName) != null : 
          "Error in module '${key}' id '${id}': required input argument '${arg.plainName}' is missing"
      }
    }

    inputs = inputs.collectEntries { name, value ->
      def par = config.allArguments.find { it.plainName == name && (it.direction == "input" || it.type == "file") }
      assert par != null : "Error in module '${key}' id '${id}': '${name}' is not a valid input argument"

      value = _checkArgumentType("input", par, value, "in module '$key' id '$id'")

      [ name, value ]
    }
  }
  return inputs
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/arguments/_processOutputValues.nf'
Map _processOutputValues(Map outputs, Map config, String id, String key) {
  if (!workflow.stubRun) {
    config.allArguments.each { arg ->
      if (arg.direction == "output" && arg.required) {
        assert outputs.containsKey(arg.plainName) && outputs.get(arg.plainName) != null : 
          "Error in module '${key}' id '${id}': required output argument '${arg.plainName}' is missing"
      }
    }

    outputs = outputs.collectEntries { name, value ->
      def par = config.allArguments.find { it.plainName == name && it.direction == "output" }
      assert par != null : "Error in module '${key}' id '${id}': '${name}' is not a valid output argument"
      
      value = _checkArgumentType("output", par, value, "in module '$key' id '$id'")
      
      [ name, value ]
    }
  }
  return outputs
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/IDChecker.nf'
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
// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/_checkUniqueIds.nf'

/**
 * Check if the ids are unique across parameter sets
 *
 * @param parameterSets a list of parameter sets.
 */
private void _checkUniqueIds(List<Tuple2<String, Map<String, Object>>> parameterSets) {
  def ppIds = parameterSets.collect{it[0]}
  assert ppIds.size() == ppIds.unique().size() : "All argument sets should have unique ids. Detected ids: $ppIds"
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/_getChild.nf'

// helper functions for reading params from file //
def _getChild(parent, child) {
  if (child.contains("://") || java.nio.file.Paths.get(child).isAbsolute()) {
    child
  } else {
    def parentAbsolute = java.nio.file.Paths.get(parent).toAbsolutePath().toString()
    parentAbsolute.replaceAll('/[^/]*$', "/") + child
  }
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/_parseParamList.nf'
/**
  * Figure out the param list format based on the file extension
  *
  * @param param_list A String containing the path to the parameter list file.
  *
  * @return A String containing the format of the parameter list file.
  */
def _paramListGuessFormat(param_list) {
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


/**
  * Read the param list
  * 
  * @param param_list One of the following:
  *   - A String containing the path to the parameter list file (csv, json or yaml),
  *   - A yaml blob of a list of maps (yaml_blob),
  *   - Or a groovy list of maps (asis).
  * @param config A Map of the Viash configuration.
  * 
  * @return A List of Maps containing the parameters.
  */
def _parseParamList(param_list, Map config) {
  // first determine format by extension
  def paramListFormat = _paramListGuessFormat(param_list)

  def paramListPath = (paramListFormat != "asis" && paramListFormat != "yaml_blob") ?
    file(param_list, hidden: true) :
    null

  // get the correct parser function for the detected params_list format
  def paramSets = []
  if (paramListFormat == "asis") {
    paramSets = param_list
  } else if (paramListFormat == "yaml_blob") {
    paramSets = readYamlBlob(param_list)
  } else if (paramListFormat == "yaml") {
    paramSets = readYaml(paramListPath)
  } else if (paramListFormat == "json") {
    paramSets = readJson(paramListPath)
  } else if (paramListFormat == "csv") {
    paramSets = readCsv(paramListPath)
  } else {
    error "Format of provided --param_list not recognised.\n" +
    "Found: '$paramListFormat'.\n" +
    "Expected: a csv file, a json file, a yaml file,\n" +
    "a yaml blob or a groovy list of maps."
  }

  // data checks
  assert paramSets instanceof List: "--param_list should contain a list of maps"
  for (value in paramSets) {
    assert value instanceof Map: "--param_list should contain a list of maps"
  }

  // id is argument
  def idIsArgument = config.allArguments.any{it.plainName == "id"}

  // Reformat from List<Map> to List<Tuple2<String, Map>> by adding the ID as first element of a Tuple2
  paramSets = paramSets.collect({ data ->
    def id = data.id
    if (!idIsArgument) {
      data = data.findAll{k, v -> k != "id"}
    }
    [id, data]
  })

  // Split parameters with 'multiple: true'
  paramSets = paramSets.collect({ id, data ->
    data = _splitParams(data, config)
    [id, data]
  })
  
  // The paths of input files inside a param_list file may have been specified relatively to the
  // location of the param_list file. These paths must be made absolute.
  if (paramListPath) {
    paramSets = paramSets.collect({ id, data ->
      def new_data = data.collectEntries{ parName, parValue ->
        def par = config.allArguments.find{it.plainName == parName}
        if (par && par.type == "file" && par.direction == "input") {
          if (parValue instanceof Collection) {
            parValue = parValue.collectMany{path -> 
              def x = _resolveSiblingIfNotAbsolute(path, paramListPath)
              x instanceof Collection ? x : [x]
            }
          } else {
            parValue = _resolveSiblingIfNotAbsolute(parValue, paramListPath) 
          }
        }
        [parName, parValue]
      }
      [id, new_data]
    })
  }

  return paramSets
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/_splitParams.nf'
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
    def parameterSettings = config.allArguments.find({it.plainName == parName})

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

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/channelFromParams.nf'
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
  // todo: fetch key from run args
  def key_ = config.name
  
  /* parse regular parameters (not in param_list)  */
  /*************************************************/
  def globalParams = config.allArguments
    .findAll { params.containsKey(it.plainName) }
    .collectEntries { [ it.plainName, params[it.plainName] ] }
  def globalID = params.get("id", null)

  /* process params_list arguments */
  /*********************************/
  def paramList = params.containsKey("param_list") && params.param_list != null ?
    params.param_list : []
  // if (paramList instanceof String) {
  //   paramList = [paramList]
  // }
  // def paramSets = paramList.collectMany{ _parseParamList(it, config) }
  // TODO: be able to process param_list when it is a list of strings
  def paramSets = _parseParamList(paramList, config)
  if (paramSets.isEmpty()) {
    paramSets = [[null, [:]]]
  }

  /* combine arguments into channel */
  /**********************************/
  def processedParams = paramSets.indexed().collect{ index, tup ->
    // Process ID
    def id = tup[0] ?: globalID
  
    if (workflow.stubRun && !id) {
      // if stub run, explicitly add an id if missing
      id = "stub${index}"
    }
    assert id != null: "Each parameter set should have at least an 'id'"

    // Process params
    def parValues = globalParams + tup[1]
    // // Remove parameters which are null, if the default is also null
    // parValues = parValues.collectEntries{paramName, paramValue ->
    //   parameterSettings = config.functionality.allArguments.find({it.plainName == paramName})
    //   if ( paramValue != null || parameterSettings.get("default", null) != null ) {
    //     [paramName, paramValue]
    //   }
    // }
    parValues = parValues.collectEntries { name, value ->
      def par = config.allArguments.find { it.plainName == name && (it.direction == "input" || it.type == "file") }
      assert par != null : "Error in module '${key_}' id '${id}': '${name}' is not a valid input argument"

      if (par == null) {
        return [:]
      }
      value = _checkArgumentType("input", par, value, "in module '$key_' id '$id'")

      [ name, value ]
    }

    [id, parValues]
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
  def processedParams = _paramsToParamSets(params, config)
  return Channel.fromList(processedParams)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/checkUniqueIds.nf'
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
// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/preprocessInputs.nf'
// This helper file will be deprecated soon
preprocessInputsDeprecationWarningPrinted = false

def preprocessInputsDeprecationWarning() {
  if (!preprocessInputsDeprecationWarningPrinted) {
    preprocessInputsDeprecationWarningPrinted = true
    System.err.println("Warning: preprocessInputs() is deprecated and will be removed in Viash 0.9.0.")
  }
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

  def config = args.config
  assert config instanceof Map : 
    "Error in preprocessInputs: config must be a map. " +
    "Expected class: Map. Found: config.getClass() is ${config.getClass()}"
  def key_ = args.key ?: config.name

  // Get different parameter types (used throughout this function)
  def defaultArgs = config.allArguments
    .findAll { it.containsKey("default") }
    .collectEntries { [ it.plainName, it.default ] }

  map { tup ->
    def id = tup[0]
    def data = tup[1]
    def passthrough = tup.drop(2)

    def new_data = (defaultArgs + data).collectEntries { name, value ->
      def par = config.allArguments.find { it.plainName == name && (it.direction == "input" || it.type == "file") }
      
      if (par != null) {
        value = _checkArgumentType("input", par, value, "in module '$key_' id '$id'")
      }

      [ name, value ]
    }

    [ id, new_data ] + passthrough
  }
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/runComponents.nf'
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
  log.warn("runComponents is deprecated, use runEach instead")
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

      def filter_ch = filter_
        ? input_ch | filter{tup ->
          filter_(tup[0], tup[1], comp_config)
        }
        : input_ch
      def id_ch = id_
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
      def data_ch = id_ch | map{tup ->
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
      def out_ch = data_ch
        | comp_.run(
          auto: (args.auto ?: [:]) + [simplifyInput: false, simplifyOutput: false]
        )
      def post_ch = toState_
        ? out_ch | map{tup ->
          def output = tup[1]
          def old_state = tup[2]
          def new_state = null
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/runEach.nf'
/**
 * Run a list of components on a stream of data.
 * 
 * @param components: list of Viash VDSL3 modules to run
 * @param fromState: a closure, a map or a list of keys to extract from the input data.
 *   If a closure, it will be called with the id, the data and the component itself.
 * @param toState: a closure, a map or a list of keys to extract from the output data
 *   If a closure, it will be called with the id, the output data, the old state and the component itself.
 * @param filter: filter function to apply to the input.
 *   It will be called with the id, the data and the component itself.
 * @param id: id to use for the output data
 *   If a closure, it will be called with the id, the data and the component itself.
 * @param auto: auto options to pass to the components
 *
 * @return: a workflow that runs the components
 **/
def runEach(Map args) {
  assert args.components: "runEach should be passed a list of components to run"

  def components_ = args.components
  if (components_ !instanceof List) {
    components_ = [ components_ ]
  }
  assert components_.size() > 0: "pass at least one component to runEach"

  def fromState_ = args.fromState
  def toState_ = args.toState
  def filter_ = args.filter
  def runIf_ = args.runIf
  def id_ = args.id

  assert !runIf_ || runIf_ instanceof Closure: "runEach: must pass a Closure to runIf."

  workflow runEachWf {
    take: input_ch
    main:

    // generate one channel per method
    out_chs = components_.collect{ comp_ ->
      def filter_ch = filter_
        ? input_ch | filter{tup ->
          filter_(tup[0], tup[1], comp_)
        }
        : input_ch
      def id_ch = id_
        ? filter_ch | map{tup ->
          def new_id = id_
          if (new_id instanceof Closure) {
            new_id = new_id(tup[0], tup[1], comp_)
          }
          assert new_id instanceof String : "Error in runEach: id should be a String or a Closure that returns a String. Expected: id instanceof String. Found: ${new_id.getClass()}"
          [new_id] + tup.drop(1)
        }
        : filter_ch
      def chPassthrough = null
      def chRun = null
      if (runIf_) {
        def idRunIfBranch = id_ch.branch{ tup ->
          run: runIf_(tup[0], tup[1], comp_)
          passthrough: true
        }
        chPassthrough = idRunIfBranch.passthrough
        chRun = idRunIfBranch.run
      } else {
        chRun = id_ch
        chPassthrough = Channel.empty()
      }
      def data_ch = chRun | map{tup ->
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
            new_data = fromState_(tup[0], new_data, comp_)
          }
          tup.take(1) + [new_data] + tup.drop(1)
        }
      def out_ch = data_ch
        | comp_.run(
          auto: (args.auto ?: [:]) + [simplifyInput: false, simplifyOutput: false]
        )
      def post_ch = toState_
        ? out_ch | map{tup ->
          def output = tup[1]
          def old_state = tup[2]
          def new_state = null
          if (toState_ instanceof Map) {
            new_state = old_state + toState_.collectEntries{ key0, key1 ->
              [key0, output[key1]]
            }
          } else if (toState_ instanceof List) {
            new_state = old_state + toState_.collectEntries{ key ->
              [key, output[key]]
            }
          } else if (toState_ instanceof Closure) {
            new_state = toState_(tup[0], output, old_state, comp_)
          }
          [tup[0], new_state] + tup.drop(3)
        }
        : out_ch

      def return_ch = post_ch
        | concat(chPassthrough)
      
      return_ch
    }

    // mix all results
    output_ch =
      (out_chs.size == 1)
        ? out_chs[0]
        : out_chs[0].mix(*out_chs.drop(1))

    emit: output_ch
  }

  return runEachWf
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/channel/safeJoin.nf'
/**
 * Join sourceChannel to targetChannel
 * 
 * This function joins the sourceChannel to the targetChannel. 
 * However, each id in the targetChannel must be present in the
 * sourceChannel. If _meta.join_id exists in the targetChannel, that is 
 * used as an id instead. If the id doesn't match any id in the sourceChannel,
 * an error is thrown.
 */

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
      
      if (!sourceIDs.contains(id)) {
        error (
          "Error in module '${key}' when merging output with original state.\n" +
          "  Reason: output with id '${id}' could not be joined with source channel.\n" +
          "    If the IDs in the output channel differ from the input channel,\n" + 
          "    please set `tup[1]._meta.join_id to the original ID.\n" +
          "  Original IDs in input channel: ['${sourceIDs.getItems().join("', '")}'].\n" + 
          "  Unexpected ID in the output channel: '${id}'.\n" +
          "  Example input event: [\"id\", [input: file(...)]],\n" +
          "  Example output event: [\"newid\", [output: file(...), _meta: [join_id: \"id\"]]]"
        )
      }
      // TODO: add link to our documentation on how to fix this

      tup
    }
  
  sourceCheck.cross(targetChannel)
    | map{ left, right ->
      right + left.drop(1)
    }
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/config/_processArgument.nf'
def _processArgument(arg) {
  arg.multiple = arg.multiple != null ? arg.multiple : false
  arg.required = arg.required != null ? arg.required : false
  arg.direction = arg.direction != null ? arg.direction : "input"
  arg.multiple_sep = arg.multiple_sep != null ? arg.multiple_sep : ";"
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/config/addGlobalParams.nf'
def addGlobalArguments(config) {
  def localConfig = [
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
          ]
          // TODO: allow multiple: true in param_list?
          // TODO: allow to specify a --param_list_regex to filter the param_list?
          // TODO: allow to specify a --param_list_from_state to remap entries in the param_list?
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/config/generateHelp.nf'
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
  def fun = config

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

// helper file: 'src/main/resources/io/viash/runners/nextflow/config/processConfig.nf'
def processConfig(config) {
  // set defaults for arguments
  config.arguments = 
    (config.arguments ?: []).collect{_processArgument(it)}

  // set defaults for argument_group arguments
  config.argument_groups =
    (config.argument_groups ?: []).collect{grp ->
      grp.arguments = (grp.arguments ?: []).collect{_processArgument(it)}
      grp
    }

  // create combined arguments list
  config.allArguments = 
    config.arguments +
    config.argument_groups.collectMany{it.arguments}

  // add missing argument groups (based on Functionality::allArgumentGroups())
  def argGroups = config.argument_groups
  if (argGroups.any{it.name.toLowerCase() == "arguments"}) {
    argGroups = argGroups.collect{ grp ->
      if (grp.name.toLowerCase() == "arguments") {
        grp = grp + [
          arguments: grp.arguments + config.arguments
        ]
      }
      grp
    }
  } else {
    argGroups = argGroups + [
      name: "Arguments",
      arguments: config.arguments
    ]
  }
  config.allArgumentGroups = argGroups

  config
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/config/readConfig.nf'

def readConfig(file) {
  def config = readYaml(file ?: moduleDir.resolve("config.vsh.yaml"))
  processConfig(config)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/_resolveSiblingIfNotAbsolute.nf'
/**
  * Resolve a path relative to the current file.
  * 
  * @param str The path to resolve, as a String.
  * @param parentPath The path to resolve relative to, as a Path.
  *
  * @return The path that may have been resovled, as a Path.
  */
def _resolveSiblingIfNotAbsolute(str, parentPath) {
  if (str !instanceof String) {
    return str
  }
  if (!_stringIsAbsolutePath(str)) {
    return parentPath.resolveSibling(str)
  } else {
    return file(str, hidden: true)
  }
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/_stringIsAbsolutePath.nf'
/**
  * Check whether a path as a string is absolute.
  *
  * In the past, we tried using `file(., relative: true).isAbsolute()`,
  * but the 'relative' option was added in 22.10.0.
  *
  * @param path The path to check, as a String.
  *
  * @return Whether the path is absolute, as a boolean.
  */
def _stringIsAbsolutePath(path) {
  def _resolve_URL_PROTOCOL = ~/^([a-zA-Z][a-zA-Z0-9]*:)?\\/.+/

  assert path instanceof String
  return _resolve_URL_PROTOCOL.matcher(path).matches()
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/collectTraces.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/deepClone.nf'
/**
  * Performs a deep clone of the given object.
  * @param x an object
  */
def deepClone(x) {
  iterateMap(x, {it instanceof Cloneable ? it.clone() : it})
}
// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/getPublishDir.nf'
def getPublishDir() {
  return params.containsKey("publish_dir") ? params.publish_dir : 
    params.containsKey("publishDir") ? params.publishDir : 
    null
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/getRootDir.nf'

// Recurse upwards until we find a '.build.yaml' file
def _findBuildYamlFile(pathPossiblySymlink) {
  def path = pathPossiblySymlink.toRealPath()
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
  def dir = _findBuildYamlFile(meta.resources_dir)
  assert dir != null: "Could not find .build.yaml in the folder structure"
  dir.getParent()
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/iterateMap.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/functions/niceView.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readCsv.nf'

def readCsv(file_path) {
  def output = []
  def inputFile = file_path !instanceof Path ? file(file_path, hidden: true) : file_path

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
        def m = removeQuote.matcher(field)
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readJson.nf'
def readJson(file_path) {
  def inputFile = file_path !instanceof Path ? file(file_path, hidden: true) : file_path
  def jsonSlurper = new groovy.json.JsonSlurper()
  jsonSlurper.parse(inputFile)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readJsonBlob.nf'
def readJsonBlob(str) {
  def jsonSlurper = new groovy.json.JsonSlurper()
  jsonSlurper.parseText(str)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readTaggedYaml.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readYaml.nf'
def readYaml(file_path) {
  def inputFile = file_path !instanceof Path ? file(file_path, hidden: true) : file_path
  def yamlSlurper = new org.yaml.snakeyaml.Yaml()
  yamlSlurper.load(inputFile)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/readYamlBlob.nf'
def readYamlBlob(str) {
  def yamlSlurper = new org.yaml.snakeyaml.Yaml()
  yamlSlurper.load(str)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/toJsonBlob.nf'
String toJsonBlob(data) {
  return groovy.json.JsonOutput.toJson(data)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/toTaggedYamlBlob.nf'
// Custom representer to modify how certain objects are represented in YAML
class CustomRepresenter extends org.yaml.snakeyaml.representer.Representer {
  Path relativizer

  class RepresentPath implements org.yaml.snakeyaml.representer.Represent {
    public String getFileName(Object obj) {
      if (obj instanceof File) {
        obj = ((File) obj).toPath();
      }
      if (obj !instanceof Path) {
        throw new IllegalArgumentException("Object: " + obj + " is not a Path or File");
      }
      def path = (Path) obj;

      if (relativizer != null) {
        return relativizer.relativize(path).toString()
      } else {
        return path.toString()
      }
    }

    public org.yaml.snakeyaml.nodes.Node representData(Object data) {
      String filename = getFileName(data);
      def tag = new org.yaml.snakeyaml.nodes.Tag("!file");
      return representScalar(tag, filename);
    }
  }
  CustomRepresenter(org.yaml.snakeyaml.DumperOptions options, Path relativizer) {
    super(options)
    this.relativizer = relativizer
    this.representers.put(sun.nio.fs.UnixPath, new RepresentPath())
    this.representers.put(Path, new RepresentPath())
    this.representers.put(File, new RepresentPath())
  }
}

String toTaggedYamlBlob(data) {
  return toRelativeTaggedYamlBlob(data, null)
}
String toRelativeTaggedYamlBlob(data, Path relativizer) {
  def options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  def representer = new CustomRepresenter(options, relativizer)
  def yaml = new org.yaml.snakeyaml.Yaml(representer, options)
  return yaml.dump(data)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/toYamlBlob.nf'
String toYamlBlob(data) {
  def options = new org.yaml.snakeyaml.DumperOptions()
  options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  options.setPrettyFlow(true)
  def yaml = new org.yaml.snakeyaml.Yaml(options)
  def cleanData = iterateMap(data, { it instanceof Path ? it.toString() : it })
  return yaml.dump(cleanData)
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/writeJson.nf'
void writeJson(data, file) {
  assert data: "writeJson: data should not be null"
  assert file: "writeJson: file should not be null"
  file.write(toJsonBlob(data))
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/readwrite/writeYaml.nf'
void writeYaml(data, file) {
  assert data: "writeYaml: data should not be null"
  assert file: "writeYaml: file should not be null"
  file.write(toYamlBlob(data))
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/states/findStates.nf'
def findStates(Map params, Map config) {
  def auto_config = deepClone(config)
  def auto_params = deepClone(params)

  auto_config = auto_config.clone()
  // override arguments
  auto_config.argument_groups = []
  auto_config.arguments = [
    [
      type: "string",
      name: "--id",
      description: "A dummy identifier",
      required: false
    ],
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
      multiple_sep: ";"
    ],
    [
      type: "string",
      name: "--settings",
      example: '{"output_dataset": "dataset.h5ad", "k": 10}',
      description: "Global arguments as a JSON glob to be passed to all components.",
      required: false
    ]
  ]
  if (!(auto_params.containsKey("id"))) {
    auto_params["id"] = "auto"
  }

  // run auto config through processConfig once more
  auto_config = processConfig(auto_config)

  workflow findStatesWf {
    helpMessage(auto_config)

    output_ch = 
      channelFromParams(auto_params, auto_config)
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
              assert split.size() == 2: "Argument 'rename_keys' should be of the form 'newKey:oldKey', or 'newKey:oldKey;newKey:oldKey' in case of multiple values"
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/states/joinStates.nf'
def joinStates(Closure apply_) {
  workflow joinStatesWf {
    take: input_ch
    main:
    output_ch = input_ch
      | toSortedList
      | filter{ it.size() > 0 }
      | map{ tups ->
        def ids = tups.collect{it[0]}
        def states = tups.collect{it[1]}
        apply_(ids, states)
      }

    emit: output_ch
  }
  return joinStatesWf
}
// helper file: 'src/main/resources/io/viash/runners/nextflow/states/publishStates.nf'
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
    def ext = path.getFileName().toString().find("\\.[^\\.]+\$") ?: ""
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
  def yamlTemplate_ = args.get("output_state", args.get("outputState", '$id.$key.state.yaml'))

  assert key_ != null : "publishStates: key must be specified"
  
  workflow publishStatesWf {
    take: input_ch
    main:
      input_ch
        | map { tup ->
          def id_ = tup[0]
          def state_ = tup[1]

          // the input files and the target output filenames
          def inputoutputFilenames_ = collectInputOutputPaths(state_, id_ + "." + key_).transpose()
          def inputFiles_ = inputoutputFilenames_[0]
          def outputFilenames_ = inputoutputFilenames_[1]

          def yamlFilename = yamlTemplate_
            .replaceAll('\\$id', id_)
            .replaceAll('\\$\\{id\\}', id_)
            .replaceAll('\\$key', key_)
            .replaceAll('\\$\\{key\\}', key_)

            // TODO: do the pathnames in state_ match up with the outputFilenames_?

          // convert state to yaml blob
          def yamlBlob_ = toRelativeTaggedYamlBlob([id: id_] + state_, java.nio.file.Paths.get(yamlFilename))

          [id_, yamlBlob_, yamlFilename, inputFiles_, outputFilenames_]
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
        [
          "[ -d \"\$(dirname '${outfile.toString()}')\" ] || mkdir -p \"\$(dirname '${outfile.toString()}')\"",
          "cp -r '${infile.toString()}' '${outfile.toString()}'"
        ]
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
  def config = args.get("config")
  assert config != null : "publishStatesByConfig: config must be specified"

  def key_ = args.get("key", config.name)
  assert key_ != null : "publishStatesByConfig: key must be specified"
  
  workflow publishStatesSimpleWf {
    take: input_ch
    main:
      input_ch
        | map { tup ->
          def id_ = tup[0]
          def state_ = tup[1] // e.g. [output: new File("myoutput.h5ad"), k: 10]
          def origState_ = tup[2] // e.g. [output: '$id.$key.foo.h5ad']

          // TODO: allow overriding the state.yaml template
          // TODO TODO: if auto.publish == "state", add output_state as an argument
          def yamlTemplate = params.containsKey("output_state") ? params.output_state : '$id.$key.state.yaml'
          def yamlFilename = yamlTemplate
            .replaceAll('\\$id', id_)
            .replaceAll('\\$\\{id\\}', id_)
            .replaceAll('\\$key', key_)
            .replaceAll('\\$\\{key\\}', key_)
          def yamlDir = java.nio.file.Paths.get(yamlFilename).getParent()

          // the processed state is a list of [key, value, inputPath, outputFilename] tuples, where
          //   - key is a String
          //   - value is any object that can be serialized to a Yaml (so a String/Integer/Long/Double/Boolean, a List, a Map, or a Path)
          //   - inputPath is a List[Path]
          //   - outputFilename is a List[String]
          //   - (key, value) are the tuples that will be saved to the state.yaml file
          //   - (inputPath, outputFilename) are the files that will be copied from src to dest (relative to the state.yaml)
          def processedState =
            config.allArguments
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
                  return [[key: plainName_, value: value, inputPath: [], outputFilename: []]]
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
                def filename = filenameTemplate
                  .replaceAll('\\$id', id_)
                  .replaceAll('\\$\\{id\\}', id_)
                  .replaceAll('\\$key', key_)
                  .replaceAll('\\$\\{key\\}', key_)
                if (par.multiple) {
                  // if the parameter is multiple: true, the filename
                  // should contain a wildcard '*' that is replaced with
                  // the index of the file
                  assert filename.contains("*") : "Module '${key_}' id '${id_}': Multiple output files specified, but no wildcard '*' in the filename: ${filename}"
                  def outputPerFile = value.withIndex().collect{ val, ix ->
                    def filename_ix = filename.replace("*", ix.toString())
                    def value_ = java.nio.file.Paths.get(filename_ix)
                    // if id contains a slash
                    if (yamlDir != null) {
                      value_ = yamlDir.relativize(value_)
                    }
                    def inputPath = val instanceof File ? val.toPath() : val
                    [value: value_, inputPath: inputPath, outputFilename: filename_ix]
                  }
                  def transposedOutputs = ["value", "inputPath", "outputFilename"].collectEntries{ key -> 
                    [key, outputPerFile.collect{dic -> dic[key]}]
                  }
                  return [[key: plainName_] + transposedOutputs]
                } else {
                  def value_ = java.nio.file.Paths.get(filename)
                  // if id contains a slash
                  if (yamlDir != null) {
                    value_ = yamlDir.relativize(value_)
                  }
                  def inputPath = value instanceof File ? value.toPath() : value
                  return [[key: plainName_, value: value_, inputPath: [inputPath], outputFilename: [filename]]]
                }
              }
          
          def updatedState_ = processedState.collectEntries{[it.key, it.value]}
          def inputPaths = processedState.collectMany{it.inputPath}
          def outputFilenames = processedState.collectMany{it.outputFilename}
          
          // convert state to yaml blob
          def yamlBlob_ = toTaggedYamlBlob([id: id_] + updatedState_)

          [id_, yamlBlob_, yamlFilename, inputPaths, outputFilenames]
        }
        | publishStatesProc
    emit: input_ch
  }
  return publishStatesSimpleWf
}

// helper file: 'src/main/resources/io/viash/runners/nextflow/states/setState.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/workflowFactory/processAuto.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/workflowFactory/processDirectives.nf'
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/workflowFactory/processWorkflowArgs.nf'
def processWorkflowArgs(Map args, Map defaultWfArgs, Map meta) {
  // override defaults with args
  def workflowArgs = defaultWfArgs + args

  // check whether 'key' exists
  assert workflowArgs.containsKey("key") : "Error in module '${meta.config.name}': key is a required argument"

  // if 'key' is a closure, apply it to the original key
  if (workflowArgs["key"] instanceof Closure) {
    workflowArgs["key"] = workflowArgs["key"](meta.config.name)
  }
  def key = workflowArgs["key"]
  assert key instanceof CharSequence : "Expected process argument 'key' to be a String. Found: class ${key.getClass()}"
  assert key ==~ /^[a-zA-Z_]\w*$/ : "Error in module '$key': Expected process argument 'key' to consist of only letters, digits or underscores. Found: ${key}"

  // check for any unexpected keys
  def expectedKeys = ["key", "directives", "auto", "map", "mapId", "mapData", "mapPassthrough", "filter", "runIf", "fromState", "toState", "args", "renameKeys", "debug"]
  def unexpectedKeys = workflowArgs.keySet() - expectedKeys
  assert unexpectedKeys.isEmpty() : "Error in module '$key': unexpected arguments to the '.run()' function: '${unexpectedKeys.join("', '")}'"

  // check whether directives exists and apply defaults
  assert workflowArgs.containsKey("directives") : "Error in module '$key': directives is a required argument"
  assert workflowArgs["directives"] instanceof Map : "Error in module '$key': Expected process argument 'directives' to be a Map. Found: class ${workflowArgs['directives'].getClass()}"
  workflowArgs["directives"] = processDirectives(defaultWfArgs.directives + workflowArgs["directives"])

  // check whether directives exists and apply defaults
  assert workflowArgs.containsKey("auto") : "Error in module '$key': auto is a required argument"
  assert workflowArgs["auto"] instanceof Map : "Error in module '$key': Expected process argument 'auto' to be a Map. Found: class ${workflowArgs['auto'].getClass()}"
  workflowArgs["auto"] = processAuto(defaultWfArgs.auto + workflowArgs["auto"])

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

  for (nam in ["map", "mapId", "mapData", "mapPassthrough", "filter", "runIf"]) {
    if (workflowArgs.containsKey(nam) && workflowArgs[nam]) {
      assert workflowArgs[nam] instanceof Closure : "Error in module '$key': Expected process argument '$nam' to be null or a Closure. Found: class ${workflowArgs[nam].getClass()}"
    }
  }

  // TODO: should functions like 'map', 'mapId', 'mapData', 'mapPassthrough' be deprecated as well?
  for (nam in ["map", "mapData", "mapPassthrough", "renameKeys"]) {
    if (workflowArgs.containsKey(nam) && workflowArgs[nam] != null) {
      log.warn "module '$key': workflow argument '$nam' is deprecated and will be removed in Viash 0.9.0. Please use 'fromState' and 'toState' instead."
    }
  }

  // check fromState
  workflowArgs["fromState"] = _processFromState(workflowArgs.get("fromState"), key, meta.config)

  // check toState
  workflowArgs["toState"] = _processToState(workflowArgs.get("toState"), key, meta.config)

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
    def requiredInputNames = meta.config.allArguments.findAll{it.required && it.direction == "Input"}.collect{it.plainName}
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
    def requiredOutputNames = config_.allArguments.findAll{it.required && it.direction == "Output"}.collect{it.plainName}
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

// helper file: 'src/main/resources/io/viash/runners/nextflow/workflowFactory/workflowFactory.nf'
def _debug(workflowArgs, debugKey) {
  if (workflowArgs.debug) {
    view { "process '${workflowArgs.key}' $debugKey tuple: $it"  }
  } else {
    map { it }
  }
}

// depends on: innerWorkflowFactory
def workflowFactory(Map args, Map defaultWfArgs, Map meta) {
  def workflowArgs = processWorkflowArgs(args, defaultWfArgs, meta)
  def key_ = workflowArgs["key"]
  
  workflow workflowInstance {
    take: input_

    main:
    def chModified = input_
      | checkUniqueIds([:])
      | _debug(workflowArgs, "input")
      | map { tuple ->
        tuple = deepClone(tuple)
        
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
        if (tuple[0] instanceof GString) {
          tuple[0] = tuple[0].toString()
        }
        assert tuple[0] instanceof CharSequence : 
          "Error in module '${key_}': first element of tuple in channel should be a String\n" +
          "  Example: [\"id\", [input: file('foo.txt'), arg: 10]].\n" +
          "  Found: ${tuple[0]}"
        
        // match file to input file
        if (workflowArgs.auto.simplifyInput && (tuple[1] instanceof Path || tuple[1] instanceof List)) {
          def inputFiles = meta.config.allArguments
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


    def chRun = null
    def chPassthrough = null
    if (workflowArgs.runIf) {
      def runIfBranch = chModified.branch{ tup ->
        run: workflowArgs.runIf(tup[0], tup[1])
        passthrough: true
      }
      chRun = runIfBranch.run
      chPassthrough = runIfBranch.passthrough
    } else {
      chRun = chModified
      chPassthrough = Channel.empty()
    }

    def chRunFiltered = workflowArgs.filter ?
      chRun | filter{workflowArgs.filter(it)} :
      chRun

    def chArgs = workflowArgs.fromState ? 
      chRunFiltered | map{
        def new_data = workflowArgs.fromState(it.take(2))
        [it[0], new_data]
      } :
      chRunFiltered | map {tup -> tup.take(2)}

    // fill in defaults
    def chArgsWithDefaults = chArgs
      | map { tuple ->
        def id_ = tuple[0]
        def data_ = tuple[1]

        // TODO: could move fromState to here

        // fetch default params from functionality
        def defaultArgs = meta.config.allArguments
          .findAll { it.containsKey("default") }
          .collectEntries { [ it.plainName, it.default ] }

        // fetch overrides in params
        def paramArgs = meta.config.allArguments
          .findAll { par ->
            def argKey = key_ + "__" + par.plainName
            params.containsKey(argKey)
          }
          .collectEntries { [ it.plainName, params[key_ + "__" + it.plainName] ] }
        
        // fetch overrides in data
        def dataArgs = meta.config.allArguments
          .findAll { data_.containsKey(it.plainName) }
          .collectEntries { [ it.plainName, data_[it.plainName] ] }
        
        // combine params
        def combinedArgs = defaultArgs + paramArgs + workflowArgs.args + dataArgs

        // remove arguments with explicit null values
        combinedArgs
          .removeAll{_, val -> val == null || val == "viash_no_value" || val == "force_null"}

        combinedArgs = _processInputValues(combinedArgs, meta.config, id_, key_)

        [id_, combinedArgs] + tuple.drop(2)
      }

    // TODO: move some of the _meta.join_id wrangling to the safeJoin() function.
    def chInitialOutput = chArgsWithDefaults
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
        def join_id = meta_.join_id ?: id_
        
        // remove metadata
        output_ = output_.findAll{k, v -> k != "_meta"}

        // check value types
        output_ = _processOutputValues(output_, meta.config, id_, key_)

        // simplify output if need be
        if (workflowArgs.auto.simplifyOutput && output_.size() == 1) {
          output_ = output_.values()[0]
        }

        [join_id, id_, output_]
      }
      // | view{"chInitialOutput: ${it.take(3)}"}

    // join the output [prev_id, new_id, output] with the previous state [prev_id, state, ...]
    def chNewState = safeJoin(chInitialOutput, chRunFiltered, key_)
      // input tuple format: [join_id, id, output, prev_state, ...]
      // output tuple format: [join_id, id, new_state, ...]
      | map{ tup ->
        def new_state = workflowArgs.toState(tup.drop(1).take(3))
        tup.take(2) + [new_state] + tup.drop(4)
      }

    if (workflowArgs.auto.publish == "state") {
      def chPublish = chNewState
        // input tuple format: [join_id, id, new_state, ...]
        // output tuple format: [join_id, id, new_state]
        | map{ tup ->
          tup.take(3)
        }

      safeJoin(chPublish, chArgsWithDefaults, key_)
        // input tuple format: [join_id, id, new_state, orig_state, ...]
        // output tuple format: [id, new_state, orig_state]
        | map { tup ->
          tup.drop(1).take(3)
      }
        | publishStatesByConfig(key: key_, config: meta.config)
    }

    // remove join_id and meta
    chReturn = chNewState
      | map { tup ->
        // input tuple format: [join_id, id, new_state, ...]
        // output tuple format: [id, new_state, ...]
        tup.drop(1)
      }
      | _debug(workflowArgs, "output")
      | concat(chPassthrough)

    emit: chReturn
  }

  def wf = workflowInstance.cloneWithName(key_)

  // add factory function
  wf.metaClass.run = { runArgs ->
    workflowFactory(runArgs, workflowArgs, meta)
  }
  // add config to module for later introspection
  wf.metaClass.config = meta.config

  return wf
}