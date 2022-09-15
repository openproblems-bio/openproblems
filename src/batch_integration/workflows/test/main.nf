nextflow.enable.dsl=2
 
targetDir = "${params.rootDir}/target/nextflow"
params.tsv = "$launchDir/src/batch_integration/datasets/params.tsv"

// import dataset loaders
include { download }       from "$targetDir/common/dataset_loader/download/main.nf"              params(params)
include { subsample }      from "$targetDir/batch_integration/datasets/subsample/main.nf"        params(params)
include { preprocessing }  from "$targetDir/batch_integration/datasets/preprocessing/main.nf"    params(params)

// import methods
include { bbknn }              from "$targetDir/batch_integration/graph/methods/bbknn/main.nf"            params(params)
include { combat }             from "$targetDir/batch_integration/graph/methods/combat/main.nf"           params(params)
include { scanorama_embed }    from "$targetDir/batch_integration/graph/methods/scanorama_embed/main.nf"  params(params)
include { scanorama_feature }    from "$targetDir/batch_integration/graph/methods/scanorama_feature/main.nf"  params(params)
include { scvi }               from "$targetDir/batch_integration/graph/methods/scvi/main.nf"             params(params)

// import metrics
include { ari }       from "$targetDir/batch_integration/graph/metrics/ari/main.nf"    params(params)
include { nmi }       from "$targetDir/batch_integration/graph/metrics/nmi/main.nf"    params(params)

/*******************************************************
*             Dataset processor workflows              *
*******************************************************/
// This workflow reads in a tsv containing some metadata about each dataset.
// For each entry in the metadata, a dataset is generated, usually by downloading
// and processing some files. The end result of each of these workflows
// should be simply a channel of [id, h5adfile, params] triplets.
//
// If the need arises, these workflows could be split off into a separate file.


workflow load_data {
    main:
        output_ = Channel.fromPath(params.tsv)
            | splitCsv(header: true, sep: "\t")
            | map { row ->
                [
                    row.name,
                    [
                        "url": row.url,
                        "name": row.name,
                        "obs_cell_type": row.label,
                        "obs_batch": row.batch,
                    ]
                ]
            }
            | download
    emit:
        output_
}

workflow process_data {

take:
      channel_in
    main:
      additional_params = Channel.fromPath(params.tsv)
        | splitCsv(header: true, sep: "\t")
        | map { [ it.name, it ] }

      subset = channel_in.join(additional_params)
        | map { id, data, additional ->
          [ id, [ input: data ] + additional ]
        }
        | subsample

      output_ = subset.join(additional_params)
        | map { id, data, additional ->
          [ id, [ input: data ] + additional ]
        }
        | preprocessing
        | join(additional_params)
        | map { id, data, additional ->
          [ id, [ input: data ] + additional + [hvg: additional.hvgs.toInteger() > 0] ]
        }

    emit:
        output_
}

/*******************************************************
*                    Main workflow                     *
*******************************************************/

workflow {
    load_data
        | process_data
        | view{ "process_data: $it" }
/**
| (bbknn & combat & scvi & scanorama_embed & scanorama_feature)
| mix
| toSortedList
        | view
        | map{ it -> [ "combined", [ input: it.collect{ it[1] } ] ] }
        | (ari & nmi)
        | extract_scores.run(
            auto: [ publish: true ]
        )
*/

}
