import numpy as np
import pandas as pd
import pybedtools
from pyensembl import EnsemblRelease
from scipy.sparse import csr_matrix
import scanpy as sc
import os
from ....tools.decorators import method


def _chrom_limit(x):
    y = x.values
    if y[-1] == "+":
        return [y[-3] - 1e5, y[-3] + 1e5]
    else:
        return [y[-2] - 1e5, y[-2] + 1e5]


def _get_annotation(adata):
    os.system("pyensembl install --release 100 --species mus_musculus")
    data = EnsemblRelease(100, species="mus_musculus")
    # get ensemble gene coordinate
    genes = []
    for i in adata.var.index.map(lambda x: x.split(".")[0]):
        try:
            gene = data.gene_by_id(i)
            genes.append(
                [
                    gene.gene_id,
                    gene.gene_name,
                    "chr%s" % gene.contig,
                    gene.start,
                    gene.end,
                    gene.strand,
                ]
            )
        except:
            genes.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    adata._var = pd.concat(
        [adata.var, pd.DataFrame(genes, index=adata.var_names)], axis=1
    )


def _atac_genes_score(adata, top_genes=500, threshold=1):
    # basic quality control
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=5)

    adata.var["mt"] = adata.var.gene_short_name.str.startswith(
        "mt-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    adata._inplace_subset_obs(adata.obs.n_genes_by_counts < 2000)
    adata._inplace_subset_obs(adata.obs.pct_counts_mt < 10)
    adata.obsm["X_tsne"] = adata.obs.loc[:, ["tsne_1", "tsne_2"]].values

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)

    # adata = adata[:, adata.var.highly_variable]
    adata._inplace_subset_var(adata.var.highly_variable)

    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    # get annotation for TSS
    _get_annotation(adata)

    # generate peak to gene weight matrix
    # remove genes without annotation
    adata._inplace_subset_var(
        (adata.var.loc[:, 2].isin(np.unique(adata.uns["mode2_var"][:, 0])))
        & (~pd.isnull(adata.var.loc[:, 2]))
    )

    # filter atac-seq matrix
    sel = np.isin(adata.uns["mode2_var"][:, 0], adata.var.loc[:, 2].unique())
    adata.uns["mode2_var"] = adata.uns["mode2_var"][sel]
    adata.obsm["mode2"] = adata.obsm["mode2"][:, sel]

    # extend tss upstream and downstream
    extend_tss = adata.var.iloc[:, -3:].apply(_chrom_limit, axis=1)

    extend_tss = pd.concat(
        [
            adata.var.loc[:, 2],
            extend_tss.map(lambda x: x[0]).astype("int32"),
            extend_tss.map(lambda x: x[1]).astype("int32"),
            pd.Series(np.arange(adata.shape[1]), index=adata.var_names),
        ],
        axis=1,
    )

    # peak summits
    peaks = pd.DataFrame(
        {
            "chr": adata.uns["mode2_var"][:, 0],
            "start": adata.uns["mode2_var"][:, 1].astype("int32"),
            "end": adata.uns["mode2_var"][:, 2].astype("int32"),
        }
    )
    summits = pd.concat(
        [
            peaks.iloc[:, 0],
            peaks.iloc[:, [1, 2]].mean(axis=1).astype("int32"),
            (peaks.iloc[:, [1, 2]].mean(axis=1) + 1).astype("int32"),
            pd.Series(np.arange(peaks.shape[0])),
        ],
        axis=1,
    )

    # overlap TSS bins with peaks
    x = pybedtools.BedTool.from_dataframe(summits)
    y = pybedtools.BedTool.from_dataframe(extend_tss)
    tss_to_peaks = x.intersect(y, wb=True, wa=True, loj=True).to_dataframe()

    # remove non-overlapped TSS and peaks
    tss_to_peaks = tss_to_peaks.loc[
        (tss_to_peaks.thickEnd != ".") | (tss_to_peaks.score != "."), :
    ]

    # weight matrix by peak to TSS distances
    tss_to_peaks["distance"] = tss_to_peaks.apply(
        lambda x: abs((int(x[5]) + int(x[6])) / 2 - int(x[1])) * 1.0 / 1e5, axis=1
    )
    tss_to_peaks["weight"] = np.exp(-0.5 - 4 * tss_to_peaks["distance"].values)

    gene_to_peak_weight = csr_matrix(
        (
            tss_to_peaks.weight.values,
            (tss_to_peaks.thickEnd.astype("int32").values, tss_to_peaks.name.values),
        ),
        shape=(adata.shape[1], adata.uns["mode2_var"].shape[0]),
    )
    adata.obsm["gene_score"] = csr_matrix.dot(
        gene_to_peak_weight, adata.obsm["mode2"].T >= threshold
    ).T


@method(
    method_name="BETA",
    paper_name="Target analysis by integration of transcriptome and ChIP-seq data with BETA",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/24263090/",
    paper_year=2013,
    code_version="1.0",
    code_url="",
)
def linear_regression_expotential_decay(adata, n_top_genes=200, threshold=1):
    _atac_genes_score(adata, top_genes=n_top_genes, threshold=threshold)
