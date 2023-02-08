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