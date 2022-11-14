nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

// // import dataset loaders
// include { download } from "$targetDir/common/dataset_loader/download/main.nf"

// // import preprocess
// include { randomize } from "$targetDir/label_projection/data_processing/randomize/main.nf"
// // for tests
// include { subsample } from "$targetDir/label_projection/data_processing/subsample/main.nf"

// // import normalization
// include { log_scran_pooling } from "$targetDir/label_projection/data_processing/normalize/scran/log_scran_pooling/main.nf"
// include { log_cpm } from "$targetDir/label_projection/data_processing/normalize/log_cpm/main.nf"

// import methods
include { all_correct } from "$targetDir/label_projection/control_methods/all_correct/main.nf"
include { majority_vote } from "$targetDir/label_projection/control_methods/majority_vote/main.nf"
include { random_labels } from "$targetDir/label_projection/control_methods/random_labels/main.nf"
include { knn_classifier } from "$targetDir/label_projection/methods/knn_classifier/main.nf"
include { mlp } from "$targetDir/label_projection/methods/mlp/main.nf"
include { logistic_regression } from "$targetDir/label_projection/methods/logistic_regression/main.nf"
// include { scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scanvi_hvg/main.nf"
// include { scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scanvi_all_genes/main.nf"
// include { scarches_scanvi_all_genes } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_all_genes/main.nf"
// include { scarches_scanvi_hvg } from "$targetDir/label_projection/methods/scvi/scarches_scanvi_hvg/main.nf"

// import metrics
include { accuracy } from "$targetDir/label_projection/metrics/accuracy/main.nf"
include { f1 } from "$targetDir/label_projection/metrics/f1/main.nf"

// import helper functions
include { extract_scores }       from "$targetDir/common/extract_scores/main.nf"


/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile] triplets.
//
// If the need arises, these workflows could be split off into a separate file.

params.tsv = "$launchDir/src/label_projection/data_processing/anndata_loader.tsv"

workflow load_data {
    main:
        output_ = Channel.fromPath(params.tsv)
            | splitCsv(header: true, sep: "\t")
            | filter{ it.name != "tabula_muris_senis_facs_lung" || it.name != "tabula_muris_senis_droplet_lung" } //TODO
            | map{ [ it.name, it ] }
            | download
    emit:
        output_
}

def mlp0 = mlp.run(
        args: [max_iter: 100, hidden_layer_sizes: 20]
)

def lr0 = logistic_regression.run(
        args: [max_iter: 100]
)

def scvi_hvg0 = scanvi_hvg.run(
    args: [n_hidden: 32, n_layers: 1, n_latent: 10, n_top_genes: 2000,
           span: 0.8, max_epochs: 1, limit_train_batches: 10, limit_val_batches: 10]
)

def scvi_allgns0 = scanvi_all_genes.run(
    args: [n_hidden: 32, n_layers: 1, n_latent: 10, n_top_genes: 2000,
           span: 0.8, max_epochs: 1, limit_train_batches: 10, limit_val_batches: 10]
)

def scarches_hvg0 = scarches_scanvi_hvg.run(
    args: [n_hidden: 32, n_layers: 1, n_latent: 10, n_top_genes: 2000,
           span: 0.8, max_epochs: 1, limit_train_batches: 10, limit_val_batches: 10]
)

def scarches_allgns0 = scarches_scanvi_all_genes.run(
    args: [n_hidden: 32, n_layers: 1, n_latent: 10, n_top_genes: 2000,
           span: 0.8, max_epochs: 1, limit_train_batches: 10, limit_val_batches: 10]
)

def f1a = f1.run(
    args: [average: "weighted"]
)

def unique_file_name(tuple) {
    return [tuple[1].baseName.replaceAll('\\.output$', ''), tuple[1]]
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    load_data
        | randomize
        | subsample.run(
            map: { [it[0], [input: it[1], even: true]] }
        )
        | (log_cpm & log_scran_pooling)
        | mix
        | map { unique_file_name(it) }
        | (knn_classifier & mlp0 & lr0 & random_labels & majority_vote & all_correct)
        | mix
        | map { unique_file_name(it) }
        | (accuracy & f1a)
        | mix
        | toSortedList
        | map{ it -> [ "combined", [ input: it.collect{ it[1] } ] ] }
        | extract_scores.run(
            args: [column_names: "dataset_id:normalization_method:method_id:metric_id:metric_value"],
            auto: [ publish: true ]
        )
}
