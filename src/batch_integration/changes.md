# Batch integration changes

Already done
* move imports to the top of the file
* remove pprint on debug
* wrote generic unit test for metrics
* merged NMI and ARI metric
* changed metric outputs to h5ad
* added .functionality.info.output_type to methods and metrics
* use anndata when scanpy is not needed
* Add 'flush=True' to all print statements
* write generic unit test for methods

TODO:
* split up api for different output types
* add file specs
* Copy descriptions from openproblems v1 for methods and metrics
* add .functionality.info to methods and metrics
* Change `split_dataset` so it uses the normalization, hvg, pca from the given object rather than recompute.
* Change `split_dataset` to output two AnnData objects: unintegrated and solution?
* Create a normalization variant:
  log_cpm -> log_cpm_batchscaled
* Rename knn_connectivities to connectivities in solution
* Simplify renaming fields in split_dataset, e.g. celltype -> label, knn_connectivities -> connectivities.

## proposed batch integration structure

* split_dataset
  - input: common dataset
  - output_unintegrated: for methods. important slots:
    .layers["normalized"]
    .obs["batch"]
    .obs["label"]? (probably not)
    .var["hvg"]
    .obsm["X_pca"]? (probably not)
  - output_solution: for metrics
    .obs["batch"]
    .obs["label"]
    .obsp["knn_connectivities"] -- could be renamed to connectivities
* methods_graph
  example: bbknn
  output_slots:
    .obsp["connectivities"]
* methods_embed
  example: scanorama_embed
  output_slots:
    .obsm["X_emb"]
* methods_feature
  example: scanorama_feature
  output_slots:
    .X
* converters
  - feature_to_embed: runs PCA on .X
  - embed_to_graph: runs sc.pp.neighbors on .obsm["X_emb"]
* metrics_graph
* metrics_embed
* metrics_feature

New todos:

- move methods and metrics to the correct folders (e.g. scanorama_embed to `methods_embed/scanorama_embed`)
- create API files for all of the above componenent types
- compare current components with v1, list the ones that are missing