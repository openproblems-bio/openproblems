from ....patch import patch_datacache

import numpy as np
import pandas as pd
import warnings


def _get_annotation(adata, retries=3):
    """Insert meta data into adata.obs."""
    from pyensembl import EnsemblRelease

    data = EnsemblRelease(
        adata.uns["release"],
        adata.uns["species"],
    )
    for _ in range(retries):
        try:
            with patch_datacache():
                data.download(overwrite=False)
                data.index(overwrite=False)
            break
        except TimeoutError:
            pass

    # get ensemble gene coordinate
    genes = []
    for i in adata.var.index.map(lambda x: x.split(".")[0]):
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    action="ignore", message="No results found for query"
                )
                gene = data.gene_by_id(i)
            genes.append(
                [
                    "chr%s" % gene.contig,
                    gene.start,
                    gene.end,
                    gene.strand,
                ]
            )
        except ValueError:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        action="ignore", message="No results found for query"
                    )
                    i = data.gene_ids_of_gene_name(i)[0]
                gene = data.gene_by_id(i)
                genes.append(
                    [
                        "chr%s" % gene.contig,
                        gene.start,
                        gene.end,
                        gene.strand,
                    ]
                )
            except (IndexError, ValueError) as e:
                print(e)
                genes.append([np.nan, np.nan, np.nan, np.nan])
    old_col = adata.var.columns.values
    adata.var = pd.concat(
        [adata.var, pd.DataFrame(genes, index=adata.var_names)], axis=1
    )
    adata.var.columns = np.hstack(
        [old_col, np.array(["chr", "start", "end", "strand"])]
    )
