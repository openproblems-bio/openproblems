#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

COMMON_DATASETS="resources/datasets/openproblems_v1"
OUTPUT_DIR="resources/label_projection/datasets/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  python << HERE
import anndata as ad
import glob
import yaml

h5ad_files = glob.glob("$COMMON_DATASETS/**.h5ad")

param_list = []

for h5ad_file in h5ad_files:
  print(f"Checking {h5ad_file}")
  adata = ad.read_h5ad(h5ad_file, backed=True)
  if "batch" in adata.obs and "celltype" in adata.obs:
    dataset_id = adata.uns["dataset_id"].replace("/", ".")
    normalization_id = adata.uns["normalization_id"]
    id = dataset_id + "." + normalization_id
    obj = {
      'id': id, 
      'input': h5ad_file,
      'dataset_id': dataset_id,
      'normalization_id': normalization_id
    }
    param_list.append(obj)

output = {
  "param_list": param_list,
  "obs_label": "celltype",
  "obs_batch": "batch",
  "seed": 123,
  "output_train": "\$id.train.h5ad",
  "output_test": "\$id.test.h5ad",
  "output_solution": "\$id.solution.h5ad"
}

with open("$params_file", "w") as file:
  yaml.dump(output, file)
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script target/nextflow/label_projection/split_dataset/main.nf \
  -profile docker \
  -resume \
  -params-file $params_file \
  --publish_dir "$OUTPUT_DIR"