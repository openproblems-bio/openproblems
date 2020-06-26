# Multimodal Data Integration

Here's a brief task description, maybe link to some seminal papers.

## API

Datasets should include matched measurements from two modalities, which are contained in `adata` and `adata.uns["mode2"]` (which itself is an `AnnData` object).

Methods should create joint matrices `adata.obsm["aligned"]` and `adata.uns["mode2"].obsm["aligned"]` which reside in a joint space.
