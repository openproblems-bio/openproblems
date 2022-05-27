nextflow.enable.dsl=2

/* For now, you need to manually specify the
 * root directory of this repository as follows.
 * (it's a nextflow limitation I'm trying to figure out 
 * how to resolve.) */
 
targetDir = "${params.rootDir}/target/nextflow"

// import dataset loaders
include { sample_dataset }       from "$targetDir/modality_alignment/datasets/sample_dataset/main.nf"     params(params)
include { scprep_csv }           from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf"         params(params)

// import methods
include { sample_method }        from "$targetDir/modality_alignment/methods/sample_method/main.nf"       params(params)
include { mnn }                  from "$targetDir/modality_alignment/methods/mnn/main.nf"                 params(params)
include { scot }                 from "$targetDir/modality_alignment/methods/scot/main.nf"                params(params)
include { harmonic_alignment }   from "$targetDir/modality_alignment/methods/harmonic_alignment/main.nf"  params(params)

// import metrics
include { knn_auc }              from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"             params(params)
include { mse }                  from "$targetDir/modality_alignment/metrics/mse/main.nf"                 params(params)

// import helper functions
include { extract_scores }       from "$targetDir/common/extract_scores/main.nf"                           params(params)


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
        output_ = Channel.fromPath("$launchDir/src/modality_alignment/datasets/datasets_scprep_csv.tsv")
            | splitCsv(header: true, sep: "\t")
            | map { row ->
                [ row.id, [ "input1": file(row.input1), "input2": file(row.input2), "id": row.id ]]
            }
            | scprep_csv
    emit:
        output_
}

workflow get_sample_datasets {
    main:
        output_ = Channel.value( [ "sample_dataset", [:] ] )
            | sample_dataset
    emit:
        output_
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    (get_sample_datasets & get_scprep_csv_datasets)
        | mix
        | (sample_method & mnn & scot & harmonic_alignment)
        | mix
        | (knn_auc & mse)
        | mix
        | toSortedList
        | map{ it -> [ "combined", [ input: it.collect{ it[1] } ] ] }
        | extract_scores.run(
            auto: [ publish: true ]
        )

}
