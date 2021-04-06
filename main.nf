nextflow.preview.dsl=2

moduleRoot="./target/nextflow/modality_alignment/"

include  { citeseq_cbmc }  from moduleRoot + 'datasets/citeseq_cbmc/main.nf'     params(params)
include  { mnn }           from moduleRoot + 'methods/mnn/main.nf'               params(params)
include  { knn_auc }       from moduleRoot + 'metrics/knn_auc/main.nf'           params(params)


workflow {
    // fetch datasets
    data_citeseq_cbmc = Channel.fromPath( "dummy" ) \
        | map{ [ "dyngen", it, params] } \
        | citeseq_cbmc

    // add more datasets here
    // data_... = ...

    // combine datasets in one channel
    datasets = data_citeseq_cbmc
    // when more datasets are available:
    // datasets = data_citeseq_cbmc.mix(data_..., data_...)

    // apply methods to datasets
    method_outputs = datasets \
        | mnn
    /*
    method_outputs = datasets \
        | (mnn & method2 & method3) \
        | mix
    */

    // apply metrics to outputs
    method_evals = method_outputs \
        | knn_auc
    /*
    method_evals = method_outputs \
        | (knn_auc & metric2 & metric3) \
        | mix
    */

    // TODO: do something with 'method_evals'
    method_evals \
        | view{ [ it[0], it[1] ] }
}
