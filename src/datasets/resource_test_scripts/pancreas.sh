#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/pancreas

set -e

mkdir -p $DATASET_DIR

# download dataset
viash run src/datasets/loaders/openproblems_v1/config.vsh.yaml -- \
    --id "pancreas" \
    --obs_celltype "celltype" \
    --obs_batch "tech" \
    --layer_counts "counts" \
    --output $DATASET_DIR/temp_dataset_full.h5ad

wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh_hm.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt
wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh_hm.txt -O $DATASET_DIR/temp_s_genes_tirosh_hm.txt
KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt $DATASET_DIR/temp_s_genes_tirosh_hm.txt | paste -sd ":" -`

# subsample
viash run src/datasets/subsample/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset_full.h5ad \
    --keep_celltype_categories "acinar:beta" \
    --keep_batch_categories "celseq:inDrop4:smarter" \
    --keep_features "$KEEP_FEATURES" \
    --output $DATASET_DIR/temp_dataset0.h5ad \
    --seed 123

# run log cpm normalisation
viash run src/datasets/normalization/log_cpm/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset0.h5ad \
    --output $DATASET_DIR/temp_dataset1.h5ad

# run pca
viash run src/datasets/processors/pca/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset1.h5ad \
    --output $DATASET_DIR/temp_dataset2.h5ad

# run log cpm normalisation
viash run src/datasets/processors/hvg/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset2.h5ad \
    --output $DATASET_DIR/dataset.h5ad

rm -r $DATASET_DIR/temp_*