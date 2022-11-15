/* usage:
| setWorkflowArguments(
  pca: [ "input": "input", "obsm_output": "obsm_pca" ]
  harmonypy: [ "obs_covariates": "obs_covariates", "obsm_input": "obsm_pca" ],
  find_neighbors: [ "obsm_input": "obsm_pca" ],
  umap: [ "output": "output" ]
)
*/

def setWorkflowArguments(Map args) {
  wfKey = args.key ?: "setWorkflowArguments"
  args.keySet().removeAll(["key"])

  
  /*
  data = [a:1, b:2, c:3]
  // args = [foo: ["a", "b"], bar: ["b"]]
  args = [foo: [a: 'a', out: "b"], bar: [in: "b"]]
  */
  
  workflow setWorkflowArgumentsInstance {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        id = tup[0]
        data = tup[1]
        passthrough = tup.drop(2)

        // determine new data
        toRemove = args.collectMany{ _, dataKeys -> 
          // dataKeys is a map but could also be a list
          dataKeys instanceof List ? dataKeys : dataKeys.values()
        }.unique()
        newData = data.findAll{!toRemove.contains(it.key)}

        // determine splitargs
        splitArgs = args.
          collectEntries{procKey, dataKeys -> 
          // dataKeys is a map but could also be a list
          newSplitData = dataKeys
            .collectEntries{ val ->
              newKey = val instanceof String ? val : val.key
              origKey = val instanceof String ? val : val.value
              [ newKey, data[origKey] ]
            }
            .findAll{it.value}
          [procKey, newSplitData]
        }

        // return output
        [ id, newData, splitArgs] + passthrough
      }

    emit:
    output_
  }

  return setWorkflowArgumentsInstance.cloneWithName(wfKey)
}

/* usage:
| getWorkflowArguments("harmonypy")
*/


def getWorkflowArguments(Map args) {
  def inputKey = args.inputKey ?: "input"
  def wfKey = "getWorkflowArguments_" + args.key
  
  workflow getWorkflowArgumentsInstance {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        id = tup[0]
        data = tup[1]
        splitArgs = tup[2].clone()
        
        passthrough = tup.drop(3)

        // try to infer arg name
        if (data !instanceof Map) {
          data = [[ inputKey, data ]].collectEntries()
        }
        newData = data + splitArgs.remove(args.key)

        [ id, newData, splitArgs] + passthrough
      }

    emit:
    output_
  }

  return getWorkflowArgumentsInstance.cloneWithName(wfKey)

}


def strictMap(Closure clos) {
  def numArgs = clos.class.methods.find{it.name == "call"}.parameterCount
  
  workflow strictMapWf {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        if (tup.size() != numArgs) {
          throw new RuntimeException("Closure does not have the same number of arguments as channel tuple.\nNumber of closure arguments: $numArgs\nChannel tuple: $tup")
        }
        clos(tup)
      }

    emit:
    output_
  }

  return strictMapWf
}

def passthroughMap(Closure clos) {
  def numArgs = clos.class.methods.find{it.name == "call"}.parameterCount
  
  workflow passthroughMapWf {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        out = clos(tup.take(numArgs))
        out + tup.drop(numArgs)
      }

    emit:
    output_
  }

  return passthroughMapWf
}

def passthroughFlatMap(Closure clos) {
  def numArgs = clos.class.methods.find{it.name == "call"}.parameterCount
  
  workflow passthroughFlatMapWf {
    take:
    input_

    main:
    output_ = input_
      | flatMap{ tup -> 
        out = clos(tup.take(numArgs))
        for (o in out) {
          o.addAll(tup.drop(numArgs))
        }
        out
      }

    emit:
    output_
  }

  return passthroughFlatMapWf
}