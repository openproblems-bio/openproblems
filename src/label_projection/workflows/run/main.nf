nextflow.enable.dsl=2

/* For now, you need to manually specify the
 * root directory of this repository as follows.
 * (it's a nextflow limitation I'm trying to figure out
 * how to resolve.) */

targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { data_loader }       from "$targetDir/common/data_loader/main.nf"     params(params)

// import preprocess
include { pancreas_preprocess }        from "$targetDir/label_projection/data/preprocess/pancreas_preprocess/main.nf"     params(params)

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
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

params.tsv = "$launchDir/src/common/data_loader/anndata_loader.tsv"

workflow load_data {
    main:
        output_ = Channel.fromPath(params.tsv)
            | splitCsv(header: true, sep: "\t")
            | map { row ->
                [ row.name, [ "url": row.url, "name": row.name ]]
            }
            | data_loader
            | pancreas_preprocess
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
