nextflow.enable.dsl=2

/* For now, you need to manually specify the
 * root directory of this repository as follows.
 * (it's a nextflow limitation I'm trying to figure out
 * how to resolve.) */

targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { data_loader }       from "$targetDir/common/data_loader/main.nf"     params(params)
// include { scprep_csv }           from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf" params(params)

// import methods
// include { sample_method }        from "$targetDir/modality_alignment/methods/sample_method/main.nf"       params(params)
// include { mnn }                  from "$targetDir/modality_alignment/methods/mnn/main.nf"                 params(params)
// include { scot }                 from "$targetDir/modality_alignment/methods/scot/main.nf"                params(params)
// include { harmonic_alignment }   from "$targetDir/modality_alignment/methods/harmonic_alignment/main.nf"  params(params)

// import metrics
// include { knn_auc }              from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"             params(params)
// include { mse }                  from "$targetDir/modality_alignment/metrics/mse/main.nf"                 params(params)

// import helper functions
// include { extract_scores }       from "$targetDir/common/extract_scores/main.nf"                           params(params)


/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

workflow load_data {
    main:
    output_ = Channel.value(params)
    | map { params -> ["new-id", ["url": params.url, "name": params.name]] }
    | data_loader
    emit:
        output_
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    load_data
    | view()

}
