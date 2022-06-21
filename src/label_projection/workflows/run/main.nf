nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { download } from "$targetDir/common/dataset_loader/download/main.nf"

// import preprocess
include { randomize } from "$targetDir/label_projection/data_processing/randomize/main.nf"

// import methods
// TODO
// import metrics
// TODO


/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

params.tsv = "$launchDir/src/common/dataset_loader/download/datasets.tsv"

workflow load_data {
    main:
        output_ = Channel.fromPath(params.tsv)
            | splitCsv(header: true, sep: "\t")
            | filter{ it.obs_celltype != "NA" && it.obs_batch != "NA" }
            | map{ [ it.name, it ] }
            | download
    emit:
        output_
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    load_data
    | randomize
    | view()
}
