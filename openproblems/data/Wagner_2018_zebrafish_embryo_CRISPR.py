from . import utils

import anndata
import scprep


@utils.loader(
    data_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112294",
    data_reference="https://doi.org/10.1126/science.aar4362",
)
def load_zebrafish_chd_tyr(test=False):
    """Download zebrafish data from GEO accession GSE112294"""

    # Information about the files to download from GEO
    sample_info = [
        ("GSM3067201", "chd", "A"),
        ("GSM3067202", "chd", "B"),
        ("GSM3067203", "chd", "C"),
        ("GSM3067204", "tyr", "A"),
        ("GSM3067205", "tyr", "B"),
        ("GSM3067206", "tyr", "C"),
    ]
    counts_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"
        "GSM3067nnn/{accession}/suppl/{accession}_{genotype}{replicate}"
        ".csv.gz"
    )
    clusters_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"
        "GSM3067nnn/{accession}/suppl/{accession}_{genotype}{replicate}_"
        "clustID.txt.gz"
    )
    cluster_names_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112294/"
        "suppl/GSE112294_ClusterNames.csv.gz"
    )

    # Download each sample and append relevant UMI / cluster info to a list
    sparse = True
    counts_matrices = []
    batch_labels = []
    metadata = []
    for accession, genotype, replicate in sample_info:
        curr_label = "{}{}".format(genotype, replicate)
        batch_labels.append(curr_label)

        # Load UMI data from GEO
        data = scprep.io.load_csv(
            counts_url.format(
                accession=accession,
                genotype=genotype,
                replicate=replicate,
            ),
            sparse=sparse,
            cell_axis="column",
        )
        counts_matrices.append(data)

        # Load cluster IDs from GEO
        clusters = scprep.io.load_csv(
            clusters_url.format(
                accession=accession, genotype=genotype, replicate=replicate
            ),
            cell_names=data.index,
            gene_names=["clusterID"],
            sparse=sparse,
        )
        metadata.append(clusters)

    # Merge data files from each sample
    data, sample_labels = scprep.utils.combine_batches(
        counts_matrices, batch_labels, append_to_cell_names=True
    )
    metadata, _ = scprep.utils.combine_batches(
        metadata, batch_labels, append_to_cell_names=True
    )
    metadata["sample"] = sample_labels
    genotype = []
    condition = []
    replicate = []
    for label in metadata["sample"]:
        if label.startswith("chd"):
            genotype.append("chd")
            condition.append("treatment")
        else:
            genotype.append("tyr")
            condition.append("control")
        replicate.append({"A": "1", "B": "2", "C": "3"}[label[-1]])

    metadata["genotype"] = genotype
    metadata["condition"] = condition
    metadata["replicate"] = replicate

    # Making cluster names more human-readable
    ClusterNamesMaps = scprep.io.load_csv(
        cluster_names_url, cell_names=False
    ).set_index("ClusterID")
    ClusterNamesMaps["ClusterName"] = ClusterNamesMaps["ClusterName"].str.slice(6)
    cluster_names = ClusterNamesMaps["ClusterName"].loc[metadata["clusterID"]]
    cluster_names.index = metadata.index
    metadata["cluster"] = cluster_names

    # Filtering rare genes and cells with high library size
    data = scprep.filter.filter_rare_genes(data)
    data, metadata = scprep.filter.filter_library_size(
        data, metadata, cutoff=15000, keep_cells="below"
    )

    # Removing cells with abnormally high expression of unannotated genes
    data, metadata = scprep.filter.filter_gene_set_expression(
        data, metadata, genes=["LOC101885394"], cutoff=164
    )

    # Convert dataframes to AnnData
    adata = anndata.AnnData(data, obs=metadata)

    # Sample uniformly from the "genotype" field to
    if test:
        adata = utils.subsample_even(adata, n_obs=500, even_obs="sample")
        utils.filter_genes_cells(adata)
        adata = adata[:, :100]
        utils.filter_genes_cells(adata)
    return adata
