nextflow.enable.dsl=2

// This workflow assumes that the directory from which the
// pipeline is launched is the root of the opsca repository.

/*******************************************************
*                 Import viash modules                 *
*******************************************************/

targetDir = "$launchDir/target/nextflow"

include { sample_dataset }       from "$targetDir/modality_alignment/datasets/sample_dataset/main.nf"     params(params)
include { scprep_csv }           from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf"         params(params)
include { sample_method }        from "$targetDir/modality_alignment/methods/sample_method/main.nf"       params(params)
include { mnn }                  from "$targetDir/modality_alignment/methods/mnn/main.nf"                 params(params)
include { scot }                 from "$targetDir/modality_alignment/methods/scot/main.nf"                params(params)
include { harmonic_alignment }   from "$targetDir/modality_alignment/methods/harmonic_alignment/main.nf"  params(params)
include { knn_auc }              from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"             params(params)
include { mse }                  from "$targetDir/modality_alignment/metrics/mse/main.nf"                 params(params)
include { extract_scores }       from "$targetDir/utils/extract_scores/main.nf"                           params(params)

include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

// Helper function for redefining the ids of elements in a channel
// based on its files.
def renameID = { [ it[1].baseName, it[1], it[2] ] }

/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

workflow get_scprep_csv_datasets {
    main:
        output_ = Channel.fromPath(file("$launchDir/src/modality_alignment/datasets/datasets_scprep_csv.tsv")) \
            | splitCsv(header: true, sep: "\t") \
            | map { row ->
                files = [ "input1": file(row.input1), "input2": file(row.input2) ]
                newParams = overrideParams(params, row.processor, "id", row.id)
                [ row.id, files, newParams, row ]
            } \
            | map{ overrideOptionValue(it, "scprep_csv", "compression", it[3].compression)} \
            | scprep_csv
    emit:
        output_
}

workflow get_sample_datasets {
    main:
        output_ = Channel.fromList( [[ "sample_dataset", [], params]] ) \
            | sample_dataset
    emit:
        output_
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    (get_sample_datasets & get_scprep_csv_datasets) \
        | mix \
        | (sample_method & mnn & scot & harmonic_alignment) \
        | mix | map(renameID) \
        | (knn_auc & mse) \
        | mix | map(renameID) \
        | toSortedList \
        | map{ it -> [ "combined", it.collect{ a -> a[1] }, params ] }
        | extract_scores

}
