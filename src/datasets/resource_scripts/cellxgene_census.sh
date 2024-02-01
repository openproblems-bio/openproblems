#!/bin/bash

# template for adding new datasets
#   - id: cellxgene_census/
#     species:
#     census_version: "2023-07-25"
#     obs_value_filter: "dataset_id == ''"
#     obs_batch:
#     dataset_name:
#     dataset_summary:
#     dataset_description:
#     dataset_url:
#     dataset_reference:
#     dataset_organism:

# not sure which dataset ids to use
#   - id: cellxgene_census/human_brain_atlas
#     species: homo_sapiens
#     census_version: "2023-07-25"
#     obs_value_filter: "dataset_id == ''" # <--- ?
#     obs_batch: donor_id
#     dataset_name:  Human Brain Atlas
#     dataset_summary: Single-Cell DNA Methylation and 3D Genome Human Brain Atlas
#     dataset_description: Delineating the gene regulatory programs underlying complex cell types is fundamental for understanding brain functions in health and disease. Here, we comprehensively examine human brain cell epigenomes by probing DNA methylation and chromatin conformation at single-cell resolution in over 500,000 cells from 46 brain regions. We identified 188 cell types and characterized their molecular signatures. Integrative analyses revealed concordant changes in DNA methylation, chromatin accessibility, chromatin organization, and gene expression across cell types, cortical areas, and basal ganglia structures. With these resources, we developed scMCodes that reliably predict brain cell types using their methylation status at select genomic sites. This multimodal epigenomic brain cell atlas provides new insights into the complexity of cell type-specific gene regulation in the adult human brain.
#     dataset_url: https://cellxgene.cziscience.com/collections/fdebfda9-bb9a-4b4b-97e5-651097ea07b0
#     dataset_reference: tian2023singlecell
#     dataset_organism: homo_sapiens

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: cellxgene_census/mouse_pancreas_atlas
    species: mus_musculus
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8'"
    obs_batch: donor_id
    dataset_name: Mouse Pancreatic Islet Atlas
    dataset_summary: Mouse pancreatic islet scRNA-seq atlas across sexes, ages, and stress conditions including diabetes
    dataset_description: To better understand pancreatic β-cell heterogeneity we generated a mouse pancreatic islet atlas capturing a wide range of biological conditions. The atlas contains scRNA-seq datasets of over 300,000 mouse pancreatic islet cells, of which more than 100,000 are β-cells, from nine datasets with 56 samples, including two previously unpublished datasets. The samples vary in sex, age (ranging from embryonic to aged), chemical stress, and disease status (including T1D NOD model development and two T2D models, mSTZ and db/db) together with different diabetes treatments. Additional information about data fields is available in anndata uns field 'field_descriptions' and on https://github.com/theislab/mm_pancreas_atlas_rep/blob/main/resources/cellxgene.md.
    dataset_url: https://cellxgene.cziscience.com/collections/296237e2-393d-4e31-b590-b03f74ac5070
    dataset_reference: hrovatin2023delineating
    dataset_organism: mus_musculus
  - id: cellxgene_census/hcla
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '066943a2-fdac-4b29-b348-40cede398e4e'"
    obs_batch: donor_id
    dataset_name: Human Lung Cell Atlas
    dataset_summary: An integrated cell atlas of the human lung in health and disease (core)
    dataset_description: The integrated Human Lung Cell Atlas (HLCA) represents the first large-scale, integrated single-cell reference atlas of the human lung. It consists of over 2 million cells from the respiratory tract of 486 individuals, and includes 49 different datasets. It is split into the HLCA core, and the extended or full HLCA. The HLCA core includes data of healthy lung tissue from 107 individuals, and includes manual cell type annotations based on consensus across 6 independent experts, as well as demographic, biological and technical metadata.
    dataset_url: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
    dataset_reference: sikkema2023integrated
    dataset_organism: homo_sapiens
  - id: cellxgene_census/tabula_sapiens
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '53d208b0-2cfd-4366-9866-c3c6114081bc'"
    obs_batch: [donor_id, assay]
    dataset_name: Tabula Sapiens
    dataset_summary: A multiple-organ, single-cell transcriptomic atlas of humans
    dataset_description: Tabula Sapiens is a benchmark, first-draft human cell atlas of nearly 500,000 cells from 24 organs of 15 normal human subjects. This work is the product of the Tabula Sapiens Consortium. Taking the organs from the same individual controls for genetic background, age, environment, and epigenetic effects and allows detailed analysis and comparison of cell types that are shared between tissues. Our work creates a detailed portrait of cell types as well as their distribution and variation in gene expression across tissues and within the endothelial, epithelial, stromal and immune compartments.
    dataset_url: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
    dataset_reference: consortium2022tabula
    dataset_organism: homo_sapiens
  - id: cellxgene_census/immune_cell_atlas
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '1b9d8702-5af8-4142-85ed-020eb06ec4f6'"
    obs_batch: donor_id
    dataset_name: Immune Cell Atlas
    dataset_summary: Cross-tissue immune cell analysis reveals tissue-specific features in humans
    dataset_description: Despite their crucial role in health and disease, our knowledge of immune cells within human tissues remains limited. We surveyed the immune compartment of 16 tissues from 12 adult donors by single-cell RNA sequencing and VDJ sequencing generating a dataset of ~360,000 cells. To systematically resolve immune cell heterogeneity across tissues, we developed CellTypist, a machine learning tool for rapid and precise cell type annotation. Using this approach, combined with detailed curation, we determined the tissue distribution of finely phenotyped immune cell types, revealing hitherto unappreciated tissue-specific features and clonal architecture of T and B cells. Our multitissue approach lays the foundation for identifying highly resolved immune cell types by leveraging a common reference dataset, tissue-integrated expression analysis, and antigen receptor sequencing.
    dataset_url: https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3
    dataset_reference: dominguez2022crosstissue
    dataset_organism: homo_sapiens
  - id: cellxgene_census/gtex_v9
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '4ed927e9-c099-49af-b8ce-a2652d069333'"
    obs_batch: donor_id
    dataset_name: GTEX v9
    dataset_summary: Single-nucleus cross-tissue molecular reference maps to decipher disease gene function
    dataset_description: Understanding the function of genes and their regulation in tissue homeostasis and disease requires knowing the cellular context in which genes are expressed in tissues across the body. Single cell genomics allows the generation of detailed cellular atlases in human tissues, but most efforts are focused on single tissue types. Here, we establish a framework for profiling multiple tissues across the human body at single-cell resolution using single nucleus RNA-Seq (snRNA-seq), and apply it to 8 diverse, archived, frozen tissue types (three donors per tissue). We apply four snRNA-seq methods to each of 25 samples from 16 donors, generating a cross-tissue atlas of 209,126 nuclei profiles, and benchmark them vs. scRNA-seq of comparable fresh tissues. We use a conditional variational autoencoder (cVAE) to integrate an atlas across tissues, donors, and laboratory methods. We highlight shared and tissue-specific features of tissue-resident immune cells, identifying tissue-restricted and non-restricted resident myeloid populations. These include a cross-tissue conserved dichotomy between LYVE1- and HLA class II-expressing macrophages, and the broad presence of LAM-like macrophages across healthy tissues that is also observed in disease. For rare, monogenic muscle diseases, we identify cell types that likely underlie the neuromuscular, metabolic, and immune components of these diseases, and biological processes involved in their pathology. For common complex diseases and traits analyzed by GWAS, we identify the cell types and gene modules that potentially underlie disease mechanisms. The experimental and analytical frameworks we describe will enable the generation of large-scale studies of how cellular and molecular processes vary across individuals and populations.
    dataset_url: https://cellxgene.cziscience.com/collections/a3ffde6c-7ad2-498a-903c-d58e732f7470
    dataset_reference: eraslan2022singlenucleus
    dataset_organism: homo_sapiens
  - id: cellxgene_census/human_retina_cell_atlas
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == 'd6505c89-c43d-4c28-8c4f-7351a5fd5528'"
    obs_batch: donor_id
    dataset_name: Human Retina Cell Atlas
    dataset_summary: Single cell atlas of the human retina
    dataset_description: As the light sensing part of the visual system, the human retina is composed of five classes of neuron, including photoreceptors, horizontal cells, amacrine, bipolar, and retinal ganglion cells. Each class of neuron can be further classified into subgroups with the abundance varying three orders of magnitude. Therefore, to capture all cell types in the retina and generate a complete single cell reference atlas, it is essential to scale up from currently published single cell profiling studies to improve the sensitivity. In addition, to gain a better understanding of gene regulation at single cell level, it is important to include sufficient scATAC-seq data in the reference. To fill the gap, we performed snRNA-seq and snATAC-seq for the retina from healthy donors. To further increase the size of the dataset, we then collected and incorporated publicly available datasets. All data underwent a unified preprocessing pipeline and data integration. Multiple integration methods were benchmarked by scIB, and scVI was chosen. To harness the power of multiomics, snATAC-seq datasets were also preprocessed, and scGlue was used to generate co-embeddings between snRNA-seq and snATAC-seq cells. To facilitate the public use of references, we employ CELLxGENE and UCSC Cell Browser for visualization. By combining previously published and newly generated datasets, a single cell atlas of the human retina that is composed of 2.5 million single cells from 48 donors has been generated. As a result, over 90 distinct cell types are identified based on the transcriptomics profile with the rarest cell type accounting for about 0.01% of the cell population. In addition, open chromatin profiling has been generated for over 400K nuclei via single nuclei ATAC-seq, allowing systematic characterization of cis-regulatory elements for individual cell type. Integrative analysis reveals intriguing differences in the transcriptome, chromatin landscape, and gene regulatory network among cell class, subgroup, and type. In addition, changes in cell proportion, gene expression and chromatin openness have been observed between different gender and over age. Accessible through interactive browsers, this study represents the most comprehensive reference cell atlas of the human retina to date. As part of the human cell atlas project, this resource lays the foundation for further research in understanding retina biology and diseases.
    dataset_url: https://cellxgene.cziscience.com/collections/4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e
    dataset_reference: li2023integrated
    dataset_organism: homo_sapiens
  - id: cellxgene_census/dkd
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id in ['ad0bf220-dd49-4b71-bb5c-576fee675d2b', 'e067e5ca-e53e-485f-aa8e-efd5435229c8']"
    obs_batch: donor_id
    dataset_name: Diabetic Kidney Disease
    dataset_summary: Multimodal single cell sequencing implicates chromatin accessibility and genetic background in diabetic kidney disease progression
    dataset_description: Multimodal single cell sequencing is a powerful tool for interrogating cell-specific changes in transcription and chromatin accessibility. We performed single nucleus RNA (snRNA-seq) and assay for transposase accessible chromatin sequencing (snATAC-seq) on human kidney cortex from donors with and without diabetic kidney disease (DKD) to identify altered signaling pathways and transcription factors associated with DKD. Both snRNA-seq and snATAC-seq had an increased proportion of VCAM1+ injured proximal tubule cells (PT_VCAM1) in DKD samples. PT_VCAM1 has a pro-inflammatory expression signature and transcription factor motif enrichment implicated NFkB signaling. We used stratified linkage disequilibrium score regression to partition heritability of kidney-function-related traits using publicly-available GWAS summary statistics. Cell-specific PT_VCAM1 peaks were enriched for heritability of chronic kidney disease (CKD), suggesting that genetic background may regulate chromatin accessibility and DKD progression. snATAC-seq found cell-specific differentially accessible regions (DAR) throughout the nephron that change accessibility in DKD and these regions were enriched for glucocorticoid receptor (GR) motifs. Changes in chromatin accessibility were associated with decreased expression of insulin receptor, increased gluconeogenesis, and decreased expression of the GR cytosolic chaperone, FKBP5, in the diabetic proximal tubule. Cleavage under targets and release using nuclease (CUT&RUN) profiling of GR binding in bulk kidney cortex and an in vitro model of the proximal tubule (RPTEC) showed that DAR co-localize with GR binding sites. CRISPRi silencing of GR response elements (GRE) in the FKBP5 gene body reduced FKBP5 expression in RPTEC, suggesting that reduced FKBP5 chromatin accessibility in DKD may alter cellular response to GR. We developed an open-source tool for single cell allele specific analysis (SALSA) to model the effect of genetic background on gene expression. Heterozygous germline single nucleotide variants (SNV) in proximal tubule ATAC peaks were associated with allele-specific chromatin accessibility and differential expression of target genes within cis-coaccessibility networks. Partitioned heritability of proximal tubule ATAC peaks with a predicted allele-specific effect was enriched for eGFR, suggesting that genetic background may modify DKD progression in a cell-specific manner.
    dataset_url: https://cellxgene.cziscience.com/collections/b3e2c6e3-9b05-4da9-8f42-da38a664b45b
    dataset_reference: wilson2022multimodal
    dataset_organism: homo_sapiens
  - id: cellxgene_census/hypomap
    species: mus_musculus
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == 'dbb4e1ed-d820-4e83-981f-88ef7eb55a35'"
    obs_batch: donor_id
    dataset_name: HypoMap
    dataset_summary: A unified single cell gene expression atlas of the murine hypothalamus
    dataset_description: The hypothalamus plays a key role in coordinating fundamental body functions. Despite recent progress in single-cell technologies, a unified catalogue and molecular characterization of the heterogeneous cell types and, specifically, neuronal subtypes in this brain region are still lacking. Here we present an integrated reference atlas “HypoMap” of the murine hypothalamus consisting of 384,925 cells, with the ability to incorporate new additional experiments. We validate HypoMap by comparing data collected from SmartSeq2 and bulk RNA sequencing of selected neuronal cell types with different degrees of cellular heterogeneity.
    dataset_url: https://cellxgene.cziscience.com/collections/d86517f0-fa7e-4266-b82e-a521350d6d36
    dataset_reference: steuernagel2022hypomap
    dataset_organism: mus_musculus

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: s3://openproblems-data/resources/datasets
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config \
  --labels cellxgene_census,dataset_loader
