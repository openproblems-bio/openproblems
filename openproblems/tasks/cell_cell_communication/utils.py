import scprep.run

# Helper function to obtain and convert a Ligand-Receptor resource
ligand_receptor_resource = scprep.run.RFunction(
    args="target_organism",
    body="""
        HUMAN <- 9606
        if (target_organism != HUMAN) {
            op_resource <- liana::generate_homologs(
                op_resource = liana::select_resource('Consensus')[[1]],
                target_organism = target_organism
            )
        } else {
            op_resource <- liana::select_resource('Consensus')[[1]]
        }

        dplyr::mutate(
            op_resource,
            ligand_genesymbol = source_genesymbol,
            receptor_genesymbol = target_genesymbol
        )
    """,
)
