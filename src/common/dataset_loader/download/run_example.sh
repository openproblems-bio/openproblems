bin/viash run src/common/dataset_loader/download/config.vsh.yaml -- \
  --output src/common/dataset_loader/download/resources/pancreas.h5ad \
  --url https://ndownloader.figshare.com/files/24539828 \
  --name pancreas \
  --obs_cell_type celltype \
  --obs_batch tech
