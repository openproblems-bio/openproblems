nextflow.enable.dsl=2

opscaRoot = "../../../"
targetDir = opscaRoot + "target/nextflow/"
taskDir = targetDir + "modality_alignment/"

println(baseDir)

include  { scprep_csv }     from '../../../target/nextflow/modality_alignment/datasets/scprep_csv/main.nf'       params(params)
// include  { scprep_csv }     from taskDir + 'datasets/scprep_csv/main.nf'       params(params)
include  { mnn }            from taskDir + 'methods/mnn/main.nf'               params(params)
include  { scot }           from taskDir + 'methods/scot/main.nf'              params(params)
include  { knn_auc }        from taskDir + 'metrics/knn_auc/main.nf'           params(params)
include  { extract_scores } from targetDir + 'utils/extract_scores/main.nf'    params(params)

// helper functions
// set id of event to basename of input file
def updateID = { [ it[1].baseName, it[1], it[2] ] }
// turn list of triplets into triplet of list
def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }
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


workflow {


    // idea: use tsv? -> https://github.com/biocorecrg/master_of_pores/blob/master/NanoMod/nanomod.nf#L80

    // fetch datasets
    data_scprep_csv = Channel.fromList( [
        [
            "CBMC_8K_13AB_10x", 
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz",
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz"
        ]
    ] ) \
        | map { [ it[0], [ "input1": file(it[1]), "input2": file(it[2]) ], params ]} \
        | map { overrideOptionValue(it, "scprep_csv", "id", it[0]) } \
        | scprep_csv

    // combine datasets in one channel
    datasets = data_scprep_csv

    // when more datasets are available, replace the code above with:
    // datasets = data_citeseq_cbmc.mix(data_2, data_3)

    datasets \
        | (mnn & scot) \
        | mix \
        | map(updateID) \
        | knn_auc \
        | map(updateID) \
        | toSortedList \
        | map( combineResults ) \
        | extract_scores


        /* When more metrics become available, replace '| knn_auc \' with the following:
        | (knn_auc & metric2 & metric3) \
        | mix \
        */
    
}
