nextflow.enable.dsl=2

println "projectDir : $projectDir"
println "launchDir : $launchDir"

targetDir = "$launchDir/target/nextflow"

include { scprep_csv; scprep_csv as notscprep_csv }     from "$targetDir/modality_alignment/datasets/scprep_csv/main.nf"       params(params)
include { mnn }            from "$targetDir/modality_alignment/methods/mnn/main.nf"               params(params)
include { scot }           from "$targetDir/modality_alignment/methods/scot/main.nf"              params(params)
include { knn_auc }        from "$targetDir/modality_alignment/metrics/knn_auc/main.nf"           params(params)
include { extract_scores } from "$targetDir/utils/extract_scores/main.nf"                         params(params)

// import helper functions
include { overrideOptionValue; overrideParams } from "$launchDir/src/utils/workflows/utils.nf"

def updateID = { [ it[1].baseName, it[1], it[2] ] }
def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }

workflow scprep_csv_wrap {
    take: 
        input_
    main:
        output_ = input_ \
            | filter{ it[3].processor == "scprep_csv" } \
            | map { overrideOptionValue(it, it[3].processor, "compression", it[3].compression) } \
            | scprep_csv
    emit:
        output_
}
workflow notscprep_csv_wrap {
    take: 
        input_
    main:
        output_ = input_ \
            | filter{ it[3].processor == "notscprep_csv" } \
            | notscprep_csv
    emit:
        output_
}

workflow {
    dataset_info = Channel.fromPath(params.datasets) \
        | splitCsv(header: true, sep: "\t") \
        | map { row -> 
            files = [ "input1": file(row.input1), "input2": file(row.input2) ]
            newParams = overrideParams(params, row.processor, "id", row.id)
            [ row.id, files, newParams, row ] 
        }

    dataset_info \
        | (scprep_csv_wrap & notscprep_csv_wrap) \
        | mix \
        | (mnn & scot) \
        | mix \
        | map(updateID) \
        | knn_auc \
        | map(updateID) \
        | toSortedList \
        | map( combineResults ) \
        | extract_scores
    
}