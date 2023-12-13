#!/bin/bash

# TODO: replace this with a run of the correct dataset loader once it is available

NEURIPS2021_URL="https://github.com/openproblems-bio/neurips2021_multimodal_viash/raw/main/resources_test/common"
DATASET_DIR="resources_test/common"

SUBDIR="$DATASET_DIR/bmmc_cite_starter"
mkdir -p "$SUBDIR"
wget "$NEURIPS2021_URL/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad" \
  -O "$SUBDIR/dataset_rna.h5ad"
wget "$NEURIPS2021_URL/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter.output_mod2.h5ad" \
  -O "$SUBDIR/dataset_adt.h5ad"

cat > "$SUBDIR/state.yaml" << HERE
id: bmmc_cite_starter
output_dataset_rna: !file dataset_rna.h5ad
output_dataset_other_mod: !file dataset_adt.h5ad
HERE

python - << HERE
import anndata as ad

rna = ad.read_h5ad("$SUBDIR/dataset_rna.h5ad")
mod2 = ad.read_h5ad("$SUBDIR/dataset_adt.h5ad")

rna.uns["dataset_id"] = "bmmc_cite_starter"
mod2.uns["dataset_id"] = "bmmc_cite_starter"
rna.uns["dataset_name"] = "BMMC Cite Starter"
mod2.uns["dataset_name"] = "BMMC Cite Starter"
rna.uns["dataset_url"] = "https://foo.bar"
mod2.uns["dataset_url"] = "https://foo.bar"
rna.uns["dataset_reference"] = "foo2001bar"
mod2.uns["dataset_reference"] = "foo2001bar"
rna.uns["dataset_summary"] = "summary"
mod2.uns["dataset_summary"] = "summary"
rna.uns["dataset_description"] = "description"
mod2.uns["dataset_description"] = "description"
rna.uns["dataset_organism"] = "homo_sapiens"
mod2.uns["dataset_organism"] = "homo_sapiens"

rna.write_h5ad("$SUBDIR/dataset_rna.h5ad")
mod2.write_h5ad("$SUBDIR/dataset_adt.h5ad")
HERE


SUBDIR="$DATASET_DIR/bmmc_multiome_starter"
mkdir -p "$SUBDIR"
wget "$NEURIPS2021_URL/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad" \
  -O "$SUBDIR/dataset_rna.h5ad"
wget "$NEURIPS2021_URL/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_mod2.h5ad" \
  -O "$SUBDIR/dataset_atac.h5ad"

cat > "$SUBDIR/state.yaml" << HERE
id: bmmc_multiome_starter
output_dataset_rna: !file dataset_rna.h5ad
output_dataset_other_mod: !file dataset_atac.h5ad
HERE


python - << HERE
import anndata as ad

rna = ad.read_h5ad("$SUBDIR/dataset_rna.h5ad")
mod2 = ad.read_h5ad("$SUBDIR/dataset_atac.h5ad")

rna.uns["dataset_id"] = "bmmc_multiome_starter"
mod2.uns["dataset_id"] = "bmmc_multipme_starter"
rna.uns["dataset_name"] = "BMMC Multiome Starter"
mod2.uns["dataset_name"] = "BMMC Multiome Starter"
rna.uns["dataset_url"] = "https://foo.bar"
mod2.uns["dataset_url"] = "https://foo.bar"
rna.uns["dataset_reference"] = "foo2001bar"
mod2.uns["dataset_reference"] = "foo2001bar"
rna.uns["dataset_summary"] = "summary"
mod2.uns["dataset_summary"] = "summary"
rna.uns["dataset_description"] = "description"
mod2.uns["dataset_description"] = "description"
rna.uns["dataset_organism"] = "homo_sapiens"
mod2.uns["dataset_organism"] = "homo_sapiens"

rna.write_h5ad("$SUBDIR/dataset_rna.h5ad")
mod2.write_h5ad("$SUBDIR/dataset_atac.h5ad")
HERE

# run task process dataset components
src/tasks/predict_modality/resources_test_scripts/bmmc_x_starter.sh