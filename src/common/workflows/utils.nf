// helper functions
// set id of event to basename of input file
//def updateID = { [ it[1].baseName, it[1], it[2] ] }
def updateID(triplet) { 
    return [ triplet[1].baseName, triplet[1], triplet[2] ]
 }


// turn list of triplets into triplet of list
//def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }
def combineResults(list) { 
    return [ "combined", list.collect{ a -> a[1] }, params ]
 }


// A functional approach to 'updating' a value for an option in the params Map.
def overrideOptionValue(triplet, _key, _option, _value) {
    mapCopy = triplet[2].toConfigObject().toMap() // As mentioned on https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/config/CascadingConfig.groovy

    return [
        triplet[0],
        triplet[1],
        triplet[2].collectEntries{ function, v1 ->
        (function == _key)
            ? [ (function) : v1.collectEntries{ k2, v2 ->
                (k2 == "arguments")
                    ? [ (k2) : v2.collectEntries{ k3, v3 ->
                        (k3 == _option)
                            ? [ (k3) : v3 + [ "value" : _value ] ]
                            : [ (k3) : v3 ]
                    } ]
                    : [ (k2) : v2 ]
            } ]
            : [ (function), v1 ]
        }
    ]
}

// A functional approach to 'updating' a value for an option in the params Map.
def overrideParams(params, _key, _option, _value) {
    mapCopy = params.toConfigObject().toMap() // As mentioned on https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/config/CascadingConfig.groovy

    return params.collectEntries{ function, v1 ->
        (function == _key)
        ? [ (function) : v1.collectEntries{ k2, v2 ->
            (k2 == "arguments")
                ? [ (k2) : v2.collectEntries{ k3, v3 ->
                    (k3 == _option)
                        ? [ (k3) : v3 + [ "value" : _value ] ]
                        : [ (k3) : v3 ]
                } ]
                : [ (k2) : v2 ]
        } ]
        : [ (function), v1 ]
    }

}