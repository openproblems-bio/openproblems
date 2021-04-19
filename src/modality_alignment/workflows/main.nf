nextflow.enable.dsl=2

println "projectDir : $projectDir"
println "launchDir : $launchDir"

targetDir = "$launchDir/target/nextflow"

include { scprep_csv }     from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf"       params(params)
include { mnn }            from "$targetDir/modality_alignment/methods/mnn/main.nf"               params(params)
include { scot }           from "$targetDir/modality_alignment/methods/scot/main.nf"              params(params)
include { knn_auc }        from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"           params(params)
include { extract_scores } from "$targetDir/utils/extract_scores/main.nf"                         params(params)

// import helper functions
include { overrideOptionValue } from "../../../src/utils/workflows/utils.nf"

def updateID = { [ it[1].baseName, it[1], it[2] ] }
def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }

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