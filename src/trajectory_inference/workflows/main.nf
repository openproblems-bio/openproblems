nextflow.enable.dsl=2

targetDir = "$launchDir/target/nextflow"

include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

include { download_datasets }       from "$targetDir/trajectory_inference/datasets/download_datasets/main.nf"     params(params)



/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

workflow get_dynverse_datasets {
    main:
        output_ = Channel.fromPath(file("$launchDir/src/trajectory_inference/datasets/datasets.tsv")) \
            | splitCsv(header: true, sep: "\t") \
            | map { row ->
                files =  file(row.links.download)
                newParams = overrideParams(params, "download_datasets", "id", row.id)
                [ row.id, files, newParams ]
            } \
            | download_datasets
    emit:
        output_
}
