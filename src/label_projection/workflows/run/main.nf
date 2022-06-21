nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { download } from "$targetDir/common/dataset_loader/download/main.nf"

// import preprocess
include { randomize } from "$targetDir/label_projection/data_processing/randomize/main.nf"

// import methods TODO[baseline/random_labels,
//                     knn_classifier/scran, mlp/log_cpm, mlp/scran, sklearn/classifier,
//                     scvi, logistic_regression/log_cpm, logistic_regression/scran]
include { majority_vote }            from "$targetDir/label_projection/methods/baseline/majority_vote/main.nf" params(params)
include { knn_classifier_log_cpm }   from "$targetDir/label_projection/methods/knn_classifier/knn_classifier_log_cpm/main.nf" params(params)
include { knn_classifier_scran }     from "$targetDir/label_projection/methods/knn_classifier/knn_classifier_scran/main.nf" params(params)//LOCALY OUT OF MEMORY
include { mlp_log_cpm }              from "$targetDir/label_projection/methods/mlp/mlp_log_cpm/main.nf" params(params)  //LOCALY OUT OF MEMORY
include { mlp_scran }                from "$targetDir/label_projection/methods/mlp/mlp_scran/main.nf" params(params)  //LOCALY OUT OF MEMORY

// import metrics TODO [f1]
include { accuracy }          from "$targetDir/label_projection/metrics/accuracy/main.nf" params(params)


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
params.tsv = "$launchDir/src/label_projection/resources/data_loader/fake_anndata_loader.tsv"

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
    | (majority_vote & knn_classifier_log_cpm & mlp_log_cpm)
    | mix
    | accuracy
    | view
}
