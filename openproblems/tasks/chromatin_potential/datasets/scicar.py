import numpy as np
import pandas as pd
import subprocess

from ....tools.decorators import dataset
from ....data.multimodal import scicar


def _get_annotation(adata):
    """Insert meta data into adata.obs."""
    from pyensembl import EnsemblRelease

    subprocess.call(
        [
            "pyensembl",
            "install",
            "--release",
            adata.uns["release"],
            "--species",
            adata.uns["species"],
        ]
    )
    data = EnsemblRelease(adata.uns["release"], species=adata.uns["species"])

    # get ensemble gene coordinate
    genes = []
    for i in adata.var.index.map(lambda x: x.split(".")[0]):
        try:
            gene = data.gene_by_id(i)
            genes.append(
                [
                    "chr%s" % gene.contig,
                    gene.start,
                    gene.end,
                    gene.strand,
                ]
            )
        except KeyError:
            genes.append([np.nan, np.nan, np.nan, np.nan])
    old_col = adata.var.columns.values
    adata.var = pd.concat(
        [adata.var, pd.DataFrame(genes, index=adata.var_names)], axis=1
    )
    adata.var.columns = np.hstack(
        [old_col, np.array(["chr", "start", "end", "strand"])]
    )


@dataset("sciCAR Mouse Kidney with cell clusters")
def scicar_mouse_kidney(test=False):
    adata = scicar.load_scicar_mouse_kidney(test=test)

    adata.uns["species"] = "mus_musculus"
    adata.uns["release"] = "100"

    # get annotation for TSS
    _get_annotation(adata)
    return adata
