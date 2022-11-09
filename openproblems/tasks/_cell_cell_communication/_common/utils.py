from typing import Union

import anndata
import collections
import numpy as np
import pandas as pd
import pathlib
import scipy.sparse
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


def map_gene_symbols(adata, map_filename: Union[str, pathlib.Path]):
    """Maps gene symbols from aliases to standard symbols

    Genes that map many-to-one are summed.
    Genes that map one-to-many are duplicated.

    Parameters
    ----------
    adata : anndata.AnnData
    map_filename : PathLike
        Path to csv containing gene symbol map with two columns, `gene` and `alias`

    Returns
    -------
    adata : anndata.AnnData
    """
    map_df = pd.read_csv(map_filename)
    var = adata.var.rename_axis("alias", axis=0)[[]]
    gene_match_idx = np.isin(var.index, map_df["gene"])
    var_gene_match, var = var.loc[gene_match_idx].copy(), var.loc[~gene_match_idx]
    alias_match_idx = np.isin(var.index, map_df["alias"])
    var_alias_match, var_no_map = (
        var.loc[alias_match_idx].copy(),
        var.loc[~alias_match_idx].copy(),
    )

    # fill 'gene' column
    var_alias_match = var_alias_match.reset_index().merge(
        map_df, on="alias", how="left"
    )
    var_gene_match["gene"] = var_gene_match.index
    var_no_map["gene"] = var_no_map.index

    var_dealiased = pd.concat(
        [var_gene_match.reset_index(), var_no_map.reset_index(), var_alias_match]
    )
    duplicate_idx = var_dealiased["gene"].duplicated(keep=False)
    var_dealiased_many_to_one, var_dealiased_one_to_any = (
        var_dealiased.loc[duplicate_idx],
        var_dealiased.loc[~duplicate_idx],
    )

    adata_one_to_any = adata[:, var_dealiased_one_to_any["alias"]]
    adata_one_to_any.var.index = var_dealiased_one_to_any["gene"]

    many_to_one_genes = var_dealiased_many_to_one["gene"].unique()
    many_to_one_X = []
    many_to_one_layers = collections.defaultdict(list)
    for gene in var_dealiased_many_to_one["gene"].unique():
        gene_aliases = var_dealiased_many_to_one.loc[
            var_dealiased_many_to_one["gene"] == gene, "alias"
        ]
        adata_gene = adata[:, gene_aliases]
        many_to_one_X.append(scipy.sparse.coo_matrix(adata_gene.X.sum(axis=1)))
        for layer_name, layer in adata_gene.layers.items():
            many_to_one_layers[layer_name].append(
                scipy.sparse.coo_matrix(adata_gene.X.sum(axis=1))
            )

    return anndata.AnnData(
        X=scipy.sparse.hstack([adata_one_to_any.X] + many_to_one_X).tocsr(),
        obs=adata.obs,
        var=pd.DataFrame(
            index=np.concatenate([adata_one_to_any.var.index, many_to_one_genes])
        ),
        layers={
            layer_name: scipy.sparse.hstack(
                [adata_one_to_any.layers[layer_name]] + many_to_one_layers[layer_name]
            ).tocsr()
            for layer_name in adata.layers
        },
        uns=adata.uns,
        obsm=adata.obsm,
    )
