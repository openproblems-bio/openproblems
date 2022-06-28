nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { download } from "$targetDir/common/dataset_loader/download/main.nf"

// import preprocess
include { randomize } from "$targetDir/label_projection/data_processing/randomize/main.nf"
// for tests
include { subsample } from "$targetDir/label_projection/data_processing/subsample/main.nf"

// import normalization
include { log_scran_pooling } from "$targetDir/label_projection/data_processing/normalize/scran/log_scran_pooling/main.nf"
include { log_cpm } from "$targetDir/label_projection/data_processing/normalize/log_cpm/main.nf"

// import methods TODO[baseline/random_labels]
include { majority_vote }            from "$targetDir/label_projection/methods/baseline/majority_vote/main.nf"
include { knn_classifier }   from "$targetDir/label_projection/methods/knn_classifier/main.nf"
include { mlp }              from "$targetDir/label_projection/methods/mlp/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
include { scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scanvi_hvg/main.nf"
include { scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scanvi_all_genes/main.nf"
include { scarches_scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_all_genes/main.nf"
include { scarches_scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_hvg/main.nf"

// import metrics TODO [f1]
include { accuracy }          from "$targetDir/label_projection/metrics/accuracy/main.nf"


/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

// params.tsv = "$launchDir/src/common/data_loader/anndata_loader.tsv"
params.tsv = "$launchDir/src/label_projection/data_processing/fake_anndata_loader.tsv" //for tests

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

def mlp0 = mlp.run(
        map: { [it[0], [input: it[1], max_iter: 100, hidden_layer_sizes: 20]] }
)

def lr0 = logistic_regression.run(
        map: { [it[0], [input: it[1], max_iter: 100]] }
)

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    load_data
    | randomize
    | subsample.run(
        map: { [it[0], [input: it[1],
                        celltype_categories: "0:3",
                        tech_categories: "0:-3:-2"]] }
    )
    | view
    | log_cpm
    | view
    | (majority_vote & knn_classifier & mlp0 & lr0)
    | mix
    | view
    | accuracy
    | view
}
