nextflow.enable.dsl=2

println "projectDir : $projectDir"
println "launchDir : $launchDir"

targetDir = "$launchDir/target/nextflow"

include { scprep_csv; scprep_csv as notscprep_csv }     from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf"       params(params)
include { mnn }            from "$targetDir/modality_alignment/methods/mnn/main.nf"               params(params)
include { scot }           from "$targetDir/modality_alignment/methods/scot/main.nf"              params(params)
include { knn_auc }        from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"           params(params)
include { mse }            from "$targetDir/modality_alignment/metrics/mse/main.nf"               params(params)
include { extract_scores } from "$targetDir/utils/extract_scores/main.nf"                         params(params)

// import helper functions
include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

def renameID = { [ it[1].baseName, it[1], it[2] ] }

workflow combineResults {
    take:
        input_
    main:
        output_ = input_ \
        | toSortedList \
        | map{ it -> [ "combined", it.collect{ a -> a[1] }, params ] }
    emit:
        output_
}

workflow get_scprep_csv_datasets {
    main:
        output_ = Channel.fromPath(file("$launchDir/src/modality_alignment/datasets/datasets_scprep_csv.tsv")) \
            | splitCsv(header: true, sep: "\t") \
            | map { row -> 
                files = [ "input1": file(row.input1), "input2": file(row.input2) ]
                newParams = overrideParams(params, row.processor, "id", row.id)
                [ row.id, files, newParams, row ] 
            } \
            | map{ overrideOptionValue(it, "scprep_csv", "compression", it[3].compression)} \
            | scprep_csv
    emit:
        output_
}

workflow {
    get_scprep_csv_datasets \
        | (mnn & scot) \
        | mix | map(renameID) \
        | (knn_auc & mse) \
        | mix | map(renameID) \
        | combineResults \
        | extract_scores
    
}