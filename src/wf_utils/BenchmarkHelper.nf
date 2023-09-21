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

def initializeTracer() {
  def traces = Collections.synchronizedList([])

  // add custom trace observer which stores traces in the traces object
  session.observers.add(new CustomTraceObserver(traces))

  traces
}

def writeJson(data, file) {
  assert data: "writeJson: data should not be null"
  assert file: "writeJson: file should not be null"
  file.write(groovy.json.JsonOutput.toJson(data))
}

import org.yaml.snakeyaml.DumperOptions
import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.nodes.Node
import org.yaml.snakeyaml.nodes.Tag
import org.yaml.snakeyaml.representer.Representer
import org.yaml.snakeyaml.representer.Represent



// Custom representer to modify how certain objects are represented in YAML
class CustomRepresenter extends Representer {
  class RepresentFile implements Represent {
    public Node representData(Object data) {
      File file = (File) data;
      String value = file.name;
      Tag tag = new Tag("!file");
      return representScalar(tag, value);
    }
  }
  CustomRepresenter(DumperOptions options) {
    super(options)
    this.representers.put(File, new RepresentFile())
  }
}

String toTaggedYamlBlob(Map data) {
  def options = new DumperOptions()
  options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
  def representer = new CustomRepresenter(options)
  def yaml = new Yaml(representer, options)
  return yaml.dump(data)
}


import org.yaml.snakeyaml.TypeDescription
import org.yaml.snakeyaml.constructor.AbstractConstruct
import org.yaml.snakeyaml.constructor.Constructor

// Custom constructor to modify how certain objects are parsed from YAML
class CustomConstructor extends Constructor {
  File root

  class ConstructFile extends AbstractConstruct {
    public Object construct(Node node) {
      String filename = (String) constructScalar(node);
      if (root != null) {
        return new File(root, filename);
      }
      return new File(filename);
    }
  }

  CustomConstructor(File root = null) {
    super()
    this.root = root
    // Handling !file tag and parse it back to a File type
    this.yamlConstructors.put(new Tag("!file"), new ConstructFile())
  }
}

def readTaggedYaml(File file) {
  Constructor constructor = new CustomConstructor(file.absoluteFile.parentFile)
  Yaml yaml = new Yaml(constructor)
  return yaml.load(file.text)
}

def getPublishDir() {
  return params.containsKey("publish_dir") ? params.publish_dir : 
    params.containsKey("publishDir") ? params.publishDir : 
    null
}

process publishStateProc {
  // todo: check publishpath?
  publishDir path: "${getPublishDir()}/${id}/", mode: "copy"
  tag "$id"
  input:
    tuple val(id), val(yamlBlob), path(inputFiles)
  output:
    tuple val(id), path{["state.yaml"] + inputFiles}
  script:
  """
  echo '${yamlBlob}' > state.yaml
  """
}

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

def convertPathsToFile(obj) {
  iterateMap(obj, {x ->
    if (x instanceof File) {
      return x
    } else if (x instanceof Path)  {
      return x.toFile()
    } else {
      return x
    }
  })
}
def convertFilesToPath(obj) {
  iterateMap(obj, {x ->
    if (x instanceof Path) {
      return x
    } else if (x instanceof File)  {
      return x.toPath()
    } else {
      return x
    }
  })
}

def publishState(Map args) {
  workflow publishStateWf {
    take: input_ch
    main:
      input_ch
        | map { tup ->
          def id = tup[0]
          def state = tup[1]
          def files = collectFiles(state)
          def convertedState = [id: id] + convertPathsToFile(state)
          def yamlBlob = toTaggedYamlBlob(convertedState)
          [id, yamlBlob, files]
        }
        | publishStateProc
    emit: input_ch
  }
  return publishStateWf
}


include { processConfig; helpMessage; channelFromParams; readYamlBlob } from "./WorkflowHelper.nf"


def findStates(Map params, Map config) {
  // TODO: do a deep clone of config
  def auto_config = config.clone()
  auto_config.functionality = auto_config.functionality.clone()
  // override arguments
  auto_config.functionality.argument_groups = []
  auto_config.functionality.arguments = [
    [
      type: "file",
      name: "--input_dir",
      example: "/path/to/input/directory",
      description: "Path to input directory containing the datasets to be integrated.",
      required: true
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
          def stateFiles = file("${args.input_dir}/**/state.yaml")

          // filter state files by regex
          if (args.filter) {
            stateFiles = stateFiles.findAll{ stateFile ->
              def stateFileStr = stateFile.toString()
              def matcher = stateFileStr =~ args.filter
              matcher.matches()}
          }

          // read in states
          def states = stateFiles.collect { stateFile ->
            def state_ = convertFilesToPath(readTaggedYaml(stateFile.toFile()))
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

def setState(fun) {
  workflow setStateWf {
    take: input_ch
    main:
      output_ch = input_ch
        | map { tup ->
          def id = tup[0]
          def state = tup[1]
          def unfilteredState = fun(id, state)
          def newState = unfilteredState.findAll{key, val -> val != null}
          [id, newState] + tup.drop(2)
        }
    emit: output_ch
  }
  return setStateWf
}