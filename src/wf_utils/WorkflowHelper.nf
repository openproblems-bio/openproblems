/////////////////////////////////////
// Viash Workflow helper functions //
/////////////////////////////////////

import java.util.regex.Pattern
import java.io.BufferedReader
import java.io.FileReader
import java.nio.file.Paths
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import org.yaml.snakeyaml.Yaml

// param helpers //
def paramExists(name) {
  return params.containsKey(name) && params[name] != ""
}

def assertParamExists(name, description) {
  if (!paramExists(name)) {
    exit 1, "ERROR: Please provide a --${name} parameter ${description}"
  }
}

// helper functions for reading params from file //
def getChild(parent, child) {
  if (child.contains("://") || Paths.get(child).isAbsolute()) {
    child
  } else {
    def parentAbsolute = Paths.get(parent).toAbsolutePath().toString()
    parentAbsolute.replaceAll('/[^/]*$', "/") + child
  }
}

def readCsv(file) {
  def output = []
  def inputFile = file !instanceof File ? new File(file) : file

  // todo: allow escaped quotes in string
  // todo: allow single quotes?
  def splitRegex = Pattern.compile(''',(?=(?:[^"]*"[^"]*")*[^"]*$)''')
  def removeQuote = Pattern.compile('''"(.*)"''')

  def br = new BufferedReader(new FileReader(inputFile))

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

def readJsonBlob(str) {
  def jsonSlurper = new JsonSlurper()
  jsonSlurper.parseText(str)
}

def readJson(file) {
  def inputFile = file !instanceof File ? new File(file) : file
  def jsonSlurper = new JsonSlurper()
  jsonSlurper.parse(inputFile)
}

def readYamlBlob(str) {
  def yamlSlurper = new Yaml()
  yamlSlurper.load(str)
}

def readYaml(file) {
  def inputFile = file !instanceof File ? new File(file) : file
  def yamlSlurper = new Yaml()
  yamlSlurper.load(inputFile)
}

// helper functions for reading a viash config in groovy //

// based on how Functionality.scala is implemented
def processArgument(arg) {
  arg.multiple = arg.multiple != null ? arg.multiple : false
  arg.required = arg.required != null ? arg.required : false
  arg.direction = arg.direction != null ? arg.direction : "input"
  arg.multiple_sep = arg.multiple_sep != null ? arg.multiple_sep : ":"
  arg.plainName = arg.name.replaceAll("^-*", "")

  if (arg.type == "file") {
    arg.must_exist = arg.must_exist != null ? arg.must_exist : true
    arg.create_parent = arg.create_parent != null ? arg.create_parent : true
  }

  if (arg.type == "file" && arg.direction == "output") {
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

// based on how Functionality.scala is implemented
def processArgumentGroup(argumentGroups, name, arguments) {
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

// based on how Functionality.scala is implemented
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
      processArgument(arg)
    }
  // set defaults for outputs
  config.functionality.outputs = 
    (config.functionality.outputs != null ? config.functionality.outputs : []).collect{arg ->
      arg.type = arg.type != null ? arg.type : "file"
      arg.direction = "output"
      processArgument(arg)
    }
  // set defaults for arguments
  config.functionality.arguments = 
    (config.functionality.arguments != null ? config.functionality.arguments : []).collect{arg ->
      processArgument(arg)
    }
  // set defaults for argument_group arguments
  config.functionality.argument_groups =
    (config.functionality.argument_groups != null ? config.functionality.argument_groups : []).collect{grp ->
      grp.arguments = (grp.arguments != null ? grp.arguments : []).collect{arg ->
        arg instanceof String ? arg.replaceAll("^-*", "") : processArgument(arg)
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
  def inputGroup = processArgumentGroup(argGroups, "Inputs", config.functionality.inputs)
  def outputGroup = processArgumentGroup(argGroups, "Outputs", config.functionality.outputs)
  def defaultGroup = processArgumentGroup(argGroups, "Arguments", config.functionality.arguments)
  def groupsFiltered = argGroups.findAll(gr -> !(["Inputs", "Outputs", "Arguments"].contains(gr.name)))
  config.functionality.allArgumentGroups = inputGroup + outputGroup + defaultGroup + groupsFiltered

  config
}

def readConfig(file) {
  def config = readYaml(file != null ? file : "$projectDir/config.vsh.yaml")
  processConfig(config)
}

// recursively merge two maps
def mergeMap(Map lhs, Map rhs) {
  return rhs.inject(lhs.clone()) { map, entry ->
    if (map[entry.key] instanceof Map && entry.value instanceof Map) {
      map[entry.key] = mergeMap(map[entry.key], entry.value)
    } else if (map[entry.key] instanceof Collection && entry.value instanceof Collection) {
      map[entry.key] += entry.value
    } else {
      map[entry.key] = entry.value
    }
    return map
  }
}

def addGlobalParams(config) {
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

  return processConfig(mergeMap(config, localConfig))
}

// helper functions for generating help // 

// based on io.viash.helpers.Format.wordWrap
def formatWordWrap(str, maxLength) {
  def words = str.split("\\s").toList()

  def word = null
  def line = ""
  def lines = []
  while(!words.isEmpty()) {
    word = words.pop()
    if (line.length() + word.length() + 1 <= maxLength) {
      line = line + " " + word
    } else {
      lines.add(line)
      line = word
    }
    if (words.isEmpty()) {
      lines.add(line)
    }
  }
  return lines
}

// based on Format.paragraphWrap
def paragraphWrap(str, maxLength) {
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

def generateArgumentHelp(param) {
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
    paragraphWrap("\n" + param.description.trim(), 80 - 8).join("\n        ")
  
  "\n    --" + param.plainName +
    namedPropsStr +
    descStr
}

// Based on Helper.generateHelp() in Helper.scala
def generateHelp(config) {
  def fun = config.functionality

  // PART 1: NAME AND VERSION
  def nameStr = fun.name + 
    (fun.version == null ? "" : " " + fun.version)

  // PART 2: DESCRIPTION
  def descrStr = fun.description == null ? 
    "" :
    "\n\n" + paragraphWrap(fun.description.trim(), 80).join("\n")

  // PART 3: Usage
  def usageStr = fun.usage == null ? 
    "" :
    "\n\nUsage:\n" + fun.usage.trim()

  // PART 4: Options
  def argGroupStrs = fun.allArgumentGroups.collect{argGroup ->
    def name = argGroup.name
    def descriptionStr = argGroup.description == null ?
      "" :
      "\n    " + paragraphWrap(argGroup.description.trim(), 80-4).join("\n    ") + "\n"
    def arguments = argGroup.arguments.collect{arg -> 
      arg instanceof String ? fun.allArguments.find{it.plainName == arg} : arg
    }.findAll{it != null}
    def argumentStrs = arguments.collect{param -> generateArgumentHelp(param)}
    
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

def helpMessage(config) {
  if (paramExists("help")) {
    def mergedConfig = addGlobalParams(config)
    def helpStr = generateHelp(mergedConfig)
    println(helpStr)
    exit 0
  }
}

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
              file(getChild(multiFile, path))
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

def paramsToChannel(params, config) {
  if (!viashChannelDeprecationWarningPrinted) {
    viashChannelDeprecationWarningPrinted = true
    System.err.println("Warning: paramsToChannel has deprecated in Viash 0.7.0. " +
                      "Please use a combination of channelFromParams and preprocessInputs.")
  }
  Channel.fromList(paramsToList(params, config))
}

def viashChannel(params, config) {
  if (!viashChannelDeprecationWarningPrinted) {
    viashChannelDeprecationWarningPrinted = true
    System.err.println("Warning: viashChannel has deprecated in Viash 0.7.0. " +
                      "Please use a combination of channelFromParams and preprocessInputs.")
  }
  paramsToChannel(params, config)
    | map{tup -> [tup.id, tup]}
}

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

/**
 * Check if the ids are unique across parameter sets
 *
 * @param parameterSets a list of parameter sets.
 */
private void _checkUniqueIds(List<Tuple2<String, Map<String, Object>>> parameterSets) {
  def ppIds = parameterSets.collect{it[0]}
  assert ppIds.size() == ppIds.unique().size() : "All argument sets should have unique ids. Detected ids: $ppIds"
}

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
          path !instanceof String ? path : file(getChild(relativeTo, path))
        })
      } else {
        parValue = parValue !instanceof String ? path : file(getChild(relativeTo, parValue))
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

  // Reformat from List<Map> to List<Tuple2<String, Map>> by adding the ID as first element of a Tuple2
  paramSets = paramSets.collect({ paramValues ->
    [paramValues.get("id", null), paramValues.findAll{it.key != 'id'}]
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
 * Cast parameters to the correct type as defined in the Viash config
 *
 * @param parValues A Map of input arguments.
 *
 * @return The input arguments that have been cast to the type from the viash config.
 */

private Map<String, Object> _castParamTypes(Map<String, Object> parValues, Map config) {
  // Cast the input to the correct type according to viash config
  def castParValues = parValues.collectEntries({ parName, parValue ->
    paramSettings = config.functionality.allArguments.find({it.plainName == parName})
    // dont parse parameters like publish_dir ( in which case paramSettings = null)
    parType = paramSettings ? paramSettings.get("type", null) : null
    if (parValue !instanceof Collection) {
      parValue = [parValue]
    }
    if (parType == "file" && ((paramSettings.direction != null ? paramSettings.direction : "input") == "input")) {
      parValue = parValue.collect{ path ->
        if (path !instanceof String) {
          path
        } else {
          file(path)
        }
      }
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
  def globalParamsValues = applyConfigToOneParameterSet(globalParams.findAll{it.key != 'id'}, config)

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
    assert id != null: "Each parameter set should have at least an ID."
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
