nextflow.enable.dsl=2

moduleRoot="./target/nextflow/modality_alignment/"

include  { citeseq_cbmc }   from moduleRoot + 'datasets/citeseq_cbmc/main.nf'     params(params)
include  { mnn }            from moduleRoot + 'methods/mnn/main.nf'               params(params)
include  { scot }           from moduleRoot + 'methods/scot/main.nf'              params(params)
include  { knn_auc }        from moduleRoot + 'metrics/knn_auc/main.nf'           params(params)
include  { extract_scores } from './target/nextflow/utils/extract_scores/main.nf' params(params)

workflow {
    // helper functions
    // set id of event to basename of input file
    def updateID = { [ it[1].baseName, it[1], it[2] ] }
    // turn list of triplets into triplet of list
    def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }

    // idea: use tsv? -> https://github.com/biocorecrg/master_of_pores/blob/master/NanoMod/nanomod.nf#L80

    // fetch datasets
    data_citeseq_cbmc = Channel.fromList( [
        [
            "citeseq_cbmc", 
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz",
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz"
        ]
    ] ) \
        | map { [ it[0], [ "input1": file(it[1]), "input2": file(it[2]) ], params ]} \
        | citeseq_cbmc

    // combine datasets in one channel
    datasets = data_citeseq_cbmc

    // when more datasets are available, replace the code above with:
    // datasets = data_citeseq_cbmc.mix(data_2, data_3)

    datasets \
        | (mnn & scot) \
        | mix \
        | map(updateID) \
        | knn_auc \
        | map(updateID) \
        | toSortedList \
        | map( combineResults ) \
        | extract_scores


        /* When more metrics become available, replace '| knn_auc \' with the following:
        | (knn_auc & metric2 & metric3) \
        | mix \
        */
    
}
