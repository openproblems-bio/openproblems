#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: zenodo_spatial/human_heart_myocardial_infarction_1_visium
    input_data: "https://zenodo.org/records/13328275/files/10X0018.h5ad?download=1"
    dataset_name: 10X Visium - Human Heart MI 1
    dataset_url: "https://www.nature.com/articles/s41586-022-05060-x"
    dataset_summary: Gene expression library of human heart using 10x Visium.
    dataset_description: "Frozen heart samples were embedded in OCT (Tissue-Tek) and cryosectioned (Thermo Cryostar). The 10-µm section was placed on the pre-chilled Optimization slides (Visium, 10X Genomics, PN-1000193) and the optimal lysis time was determined. The tissues were treated as recommended by 10X Genomics and the optimization procedure showed an optimal permeabilization time of 12 or 18 min of digestion and release of RNA from the tissue slide. Spatial gene expression slides (Visium, 10X Genomics, PN-1000187) were used for spatial transcriptomics following the Visium User Guides"
    dataset_reference: kuppe2022spatial
    dataset_organism: Homo sapiens
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: true

  - id: zenodo_spatial/human_heart_myocardial_infarction_2_visium
    input_data: "https://zenodo.org/records/13328275/files/10X009.h5ad?download=1"
    dataset_name: 10X Visium - Human Heart MI 2
    dataset_url: "https://www.nature.com/articles/s41586-022-05060-x"
    dataset_summary: Gene expression library of human heart using 10x Visium.
    dataset_description: "Frozen heart samples were embedded in OCT (Tissue-Tek) and cryosectioned (Thermo Cryostar). The 10-µm section was placed on the pre-chilled Optimization slides (Visium, 10X Genomics, PN-1000193) and the optimal lysis time was determined. The tissues were treated as recommended by 10X Genomics and the optimization procedure showed an optimal permeabilization time of 12 or 18 min of digestion and release of RNA from the tissue slide. Spatial gene expression slides (Visium, 10X Genomics, PN-1000187) were used for spatial transcriptomics following the Visium User Guides"
    dataset_reference: kuppe2022spatial
    dataset_organism: Homo sapiens
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: true

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: resources/datasets
remove_mitochondrial: true
HERE

# catt > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/mouse_e10_brain_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_brain_gene_25um_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Brain (E10)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E10 whole mouse embryo tissue (brain in early-stage organogenesis) profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_e10_eye_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_eye_and_nearby_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Eye (E10)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E10 whole mouse embryo tissue (eye in early-stage organogenesis) profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_e10_whole_body_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E10_whole_gene_best_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Whole Body (E10)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E10 whole mouse embryo tissue profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_e11_lower_body_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_E11_lower_body_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Lower Body (E11)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E11 whole mouse embryo tissue (lower body in early-stage organogenesis) profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_e11_1_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_GSM4364244_E11-FL-1L_gene_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Whole Body 1 (E11)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E11 whole mouse embryo tissue profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_e11_2_dbitseq
#     input_data: "https://zenodo.org/records/12785822/files/DBiT-seq_liu2020high_GSM4364245_E11-FL-2L_gene_data.h5ad?download=1"
#     dataset_name: DBiT-seq - Mouse Whole Body 2 (E11)
#     dataset_url: "https://www.cell.com/cell/fulltext/S0092-8674(20)31390-8"
#     dataset_summary: High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.
#     dataset_description: "Gene expression library of an E11 whole mouse embryo tissue profiled using DBiT-seq."
#     dataset_organism: Mus musculus
#     dataset_reference: liu2020high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# HERE

# cat > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/human_cortex_1_merfish
#     input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.250.expand.rep1_data.h5ad?download=1"
#     dataset_name: MERFISH - Human Cortex 1
#     dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
#     dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
#     dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 1 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (250 gene panel)."
#     dataset_organism: Homo sapiens
#     dataset_reference: fang2022conservation
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 100
#     remove_mitochondrial: false

#   - id: zenodo_spatial/human_cortex_2_merfish
#     input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep1_data.h5ad?download=1"
#     dataset_name: MERFISH - Human Cortex 2
#     dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
#     dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
#     dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 1 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
#     dataset_organism: Homo sapiens
#     dataset_reference: fang2022conservation
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: false

#   - id: zenodo_spatial/human_cortex_3_merfish
#     input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep2_data.h5ad?download=1"
#     dataset_name: MERFISH - Human Cortex 3
#     dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
#     dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
#     dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 2 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
#     dataset_organism: Homo sapiens
#     dataset_reference: fang2022conservation
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: false

#   - id: zenodo_spatial/human_cortex_4_merfish
#     input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_H18.06.006.MTG.4000.expand.rep3_data.h5ad?download=1"
#     dataset_name: MERFISH - Human Cortex 4
#     dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
#     dataset_summary: Spatially resolved profiling of human cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
#     dataset_description: "Spatially resolved profiling of human cerebral cortex (middle temopral gyrus) replicate 3 using multiplexed error-robust fluorescence in situ hybridization (MERFISH) (4000 gene panel)."
#     dataset_organism: Homo sapiens
#     dataset_reference: fang2022conservation
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: false

#   - id: zenodo_spatial/mouse_cortex_merfish
#     input_data: "https://zenodo.org/records/12785822/files/MERFISH_Fang2022Conservation_mouse1.AUD_TEA_VIS.242.unexpand_data.h5ad?download=1"
#     dataset_name: MERFISH - Mouse Cortex
#     dataset_url: "https://www.science.org/doi/10.1126/science.abm1741"
#     dataset_summary: Spatially resolved profiling of mouse cerebral cortex using multiplexed error-robust fluorescence in situ hybridization (MERFISH).
#     dataset_description: "Spatially resolved profiling of mouse cerebral cortex (visual cortex (VIS), auditory cortex (AUD) and temporal association area (TEa) unexpanded sections) using multiplexed error-robust fluorescence in situ hybridization (MERFISH)."
#     dataset_organism: Mus musculus
#     dataset_reference: fang2022conservation
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# HERE

# cat > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/mouse_organogenesis_seqfish
#     input_data: "https://zenodo.org/records/12785822/files/seqfish.h5ad?download=1"
#     dataset_name: Seqfish - Mouse Organogenesis
#     dataset_url: "https://www.nature.com/articles/s41587-021-01006-2"
#     dataset_summary: Single-cell spatial expression of mouse organogenesis.
#     dataset_description: "Sagittal sections from mouse embryo at the 8-12 ss was profiled by seqFISH."
#     dataset_organism: Mus musculus
#     dataset_reference: lohoff2021integration
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 10
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# remove_mitochondrial: true
# HERE

# cat > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/mouse_olfactory_bulb_puck_slideseqv2
#     input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_SlideSeqV2_Mouse_Olfactory_bulb_Puck_200127_15_data_whole.h5ad?download=1"
#     dataset_name: Slide-seqV2 - Mouse Olfactory Bulb Puck
#     dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
#     dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
#     dataset_description: "Gene expression library of mouse olfactory bulk puck profiled using Slide-seq V2."
#     dataset_reference: stickels2020highly
#     dataset_organism: Mus musculus
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 500
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_cortex_slideseqv2
#     input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_palla2021squidpy_Slide-seqV2_Mouse_Cortex_data_whole.h5ad?download=1"
#     dataset_name: Slide-seqV2 - Mouse Cortex
#     dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
#     dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
#     dataset_description: "Gene expression library of Mouse cortex profiled using Slide-seq V2."
#     dataset_reference: stickels2020highly
#     dataset_organism: Mus musculus
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 500
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_cerebellum_slideseqv2
#     input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_Cerebellum_SCP948_data_whole.h5ad?download=1"
#     dataset_name: Slide-seqV2 - Mouse Cerebellum
#     dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
#     dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
#     dataset_description: "Gene expression library of mouse cerebellum profiled using Slide-seq V2."
#     dataset_reference: stickels2020highly
#     dataset_organism: Mus musculus
#     spot_filter_min_genes: 100
#     gene_filter_min_spots: 500
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_hippocampus_puck_slideseqv2
#     input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_Hippocampus_Puck_200115_08_data_whole.h5ad?download=1"
#     dataset_name: Slide-seqV2 - Mouse Hippocampus Puck
#     dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
#     dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
#     dataset_description: "Gene expression library of mouse hippocampus puck profiled using Slide-seq V2."
#     dataset_reference: stickels2020highly
#     dataset_organism: Mus musculus
#     spot_filter_min_genes: 200
#     gene_filter_min_spots: 500
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_somatosensory_cortex_puck_slideseqv2
#     input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_Slide-seqV2_Mouse_SomatosensoryCortex_Puck_200306_03_data_whole.h5ad?download=1"
#     dataset_name: Slide-seqV2 - Mouse Somatosensory Cortex Puck
#     dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary"
#     dataset_summary: Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2.
#     dataset_description: "Gene expression library of mouse somatosensory cortex puck profiled using Slide-seq V2."
#     dataset_reference: stickels2020highly
#     dataset_organism: Mus musculus
#     spot_filter_min_genes: 200
#     gene_filter_min_spots: 500
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# HERE

# cat > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/mouse_brain_2d_zstep10_0_starmap
#     input_data: "https://zenodo.org/records/12785822/files/STARmap_Wang2018three_data_2D_zstep10_0_data.h5ad?download=1"
#     dataset_name: STARmap - Mouse Brain 1
#     dataset_url: "https://www.science.org/doi/10.1126/science.aat5691"
#     dataset_summary: Three-dimensional intact-tissue sequencing of single-cell transcriptional states.
#     dataset_description: "3D architecture of cell types in visual cortex volumes."
#     dataset_organism: Mus musculus
#     dataset_reference: wang2018three
#     spot_filter_min_genes: 1
#     gene_filter_min_spots: 1
#     remove_mitochondrial: true

#   - id: zenodo_spatial/mouse_brain_2d_zstep15_0_starmap
#     input_data: "https://zenodo.org/records/12785822/files/STARmap_Wang2018three_data_2D_zstep15_0_data.h5ad?download=1"
#     dataset_name: STARmap - Mouse Brain 2
#     dataset_url: "https://www.science.org/doi/10.1126/science.aat5691"
#     dataset_summary: Three-dimensional intact-tissue sequencing of single-cell transcriptional states.
#     dataset_description: "3D architecture of cell types in visual cortex volumes."
#     dataset_organism: Mus musculus
#     dataset_reference: wang2018three
#     spot_filter_min_genes: 1
#     gene_filter_min_spots: 1
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# HERE

# cat > "/tmp/params.yaml" << 'HERE'
# param_list:
#   - id: zenodo_spatial/drosophila_embryo_e5_6_stereoseq
#     input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_5.6.h5ad?download=1"
#     dataset_name: Stereo-seq - Drosophila embryo E5_6
#     dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
#     dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
#     dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
#     dataset_organism: Drosophila
#     dataset_reference: wang2022high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/drosophila_embryo_e6_3_stereoseq
#     input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_6.3.h5ad?download=1"
#     dataset_name: Stereo-seq - Drosophila embryo E6_3
#     dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
#     dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
#     dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
#     dataset_organism: Drosophila
#     dataset_reference: wang2022high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/drosophila_embryo_e7_stereoseq
#     input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_7.h5ad?download=1"
#     dataset_name: Stereo-seq - Drosophila embryo E7
#     dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
#     dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
#     dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
#     dataset_organism: Drosophila
#     dataset_reference: wang2022high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/drosophila_embryo_e9_1_stereoseq
#     input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_9.1.h5ad?download=1"
#     dataset_name: Stereo-seq - Drosophila embryo E9_1
#     dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
#     dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
#     dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
#     dataset_organism: Drosophila
#     dataset_reference: wang2022high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

#   - id: zenodo_spatial/drosophila_embryo_e10_stereoseq
#     input_data: "https://zenodo.org/records/12785822/files/Stereo-seq_wang2022high_E14-16h_a_count_normal_stereoseq_data_whole_time_point_10.5.h5ad?download=1"
#     dataset_name: Stereo-seq - Drosophila embryo E10
#     dataset_url: "https://www.sciencedirect.com/science/article/pii/S1534580722002465"
#     dataset_summary: Stereo-seq faithfully captures Drosophila spatial transcriptomes with high resolution.
#     dataset_description: "Drosophila has long been a successful model organism in multiple biomedical fields. Spatial gene expression patterns are critical for the understanding of complex pathways and interactions, whereas temporal gene expression changes are vital for studying highly dynamic physiological activities. Systematic studies in Drosophila are still impeded by the lack of spatiotemporal transcriptomic information. Here, utilizing spatial enhanced resolution omics-sequencing (Stereo-seq), we dissected the spatiotemporal transcriptomic changes of developing Drosophila with high resolution and sensitivity. (Data from an embryo collected 14-16 h after egg laying)"
#     dataset_organism: Drosophila
#     dataset_reference: wang2022high
#     spot_filter_min_genes: 10
#     gene_filter_min_spots: 50
#     remove_mitochondrial: true

# normalization_methods: [log_cp10k]
# output_dataset: '$id/dataset.h5ad'
# output_meta: '$id/dataset_metadata.yaml'
# output_state: '$id/state.yaml'
# output_raw: force_null
# output_normalized: force_null
# publish_dir: resources/datasets
# HERE

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
  --main-script target/nextflow/datasets/workflows/process_zenodo_spatial/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
