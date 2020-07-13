# Multimodal Data Integration

Here's a brief task description, maybe link to some seminal papers.

## API

Datasets should include matched measurements from two modalities, which are contained in `adata` and `adata.obsm["mode2"]`.

Methods should create joint matrices `adata.obsm["aligned"]` and `adata.obsm["mode2_aligned"]` which reside in a joint space.
