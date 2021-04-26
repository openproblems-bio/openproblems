nextflow.enable.dsl=2

// This workflow assumes that the directory from which the
// pipeline is launched is the root of the opsca repository.

/*******************************************************
*                 Import viash modules                 *
*******************************************************/

targetDir = "$launchDir/target/nextflow"

include { sample_dataset }       from "$targetDir/modality_alignment/datasets/sample_dataset/main.nf"     params(params)
include { sample_method }        from "$targetDir/modality_alignment/methods/sample_method/main.nf"       params(params)
include { knn_auc }              from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"             params(params)
include { mse }                  from "$targetDir/modality_alignment/metrics/mse/main.nf"                 params(params)
include { extract_scores }       from "$targetDir/utils/extract_scores/main.nf"                           params(params)

include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

// Helper function for redefining the ids of elements in a channel
// based on its files.
def renameID = { [ it[1].baseName, it[1], it[2] ] }


/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    Channel.fromList( [[ "sample_dataset", [], params]] ) \
        | sample_dataset \
        | sample_method \
        | map(renameID) \
        | (knn_auc & mse) \
        | mix | map(renameID) \
        | toSortedList \
        | map{ it -> [ "combined", it.collect{ a -> a[1] }, params ] }
        | extract_scores

}
