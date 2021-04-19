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

def updateID = { [ it[1].baseName, it[1], it[2] ] }
def combineResults = { it -> [ "combined", it.collect{ a -> a[1] }, params ] }

workflow multiMerge {
    take:
        input_
    main:
        output_ = input_ \
        | mix \
        | map{ [ it[1].baseName, it[1], it[2] ] }
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

workflow processDatasets {
    take: 
        dataset_info_
    main:
        // read 
        datasets = dataset_info_ \
            | splitCsv(header: true, sep: "\t") \
            | map { row -> 
                files = [ "input1": file(row.input1), "input2": file(row.input2) ]
                newParams = overrideParams(params, row.processor, "id", row.id)
                [ row.id, files, newParams, row ] 
            } \
            | branch {
                data_scprep_csv: it[3].processor == "scprep_csv"
                data_notscprep_csv: it[3].processor == "notscprep_csv"
            }

        // add extra scprep_csv parameters from tsv rows
        out_scprep_csv = datasets.data_scprep_csv \
            | map{ overrideOptionValue(it, "scprep_csv", "compression", it[3].compression)}
            | scprep_csv

        // process other data loaders
        out_notscprep_csv = datasets.data_notscprep_csv | notscprep_csv

        // combine multiple data loaders into a single channel
        output_ = out_scprep_csv.mix(out_notscprep_csv)
    emit:
        output_
}

workflow {
    get_scprep_csv_datasets \
        | (mnn & scot) \
        | mix \
        | map(updateID) \
        | (knn_auc & mse) \
        | mix \
        | map(updateID) \
        | toSortedList \
        | map(combineResults) \
        | extract_scores
    
}