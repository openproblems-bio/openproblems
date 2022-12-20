#!/bin/bash

# make sure folloewing command has been executed
# viash ns build -q 'common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUTPUT_DIR="resources_test/common/task_metadata"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

sha_file="$OUTPUT_DIR/input_git_sha.json"

cat <<EOT > $sha_file
[
  {
    "path": "tasks/denoising/README.md",
    "last_modified": "2022-09-20 14:26:51 -0400",
    "sha": "3fe9251ba906061b6769eed2ac9da0db5f8e26bb"
  },
  {
    "path": "tasks/denoising/__init__.py",
    "last_modified": "2022-09-30 14:49:17 +0200",
    "sha": "c97decf07adb2e3050561d6fa9ae46132be07bef"
  },
  {
    "path": "tasks/denoising/api.py",
    "last_modified": "2022-10-21 13:56:15 -0400",
    "sha": "b460ecb183328c857cbbf653488f522a4034a61c"
  },
  {
    "path": "tasks/denoising/datasets/__init__.py",
    "last_modified": "2022-11-23 10:32:02 -0500",
    "sha": "725ff0c46140aaa6bbded68646256f64bc63df6d"
  },
  {
    "path": "tasks/denoising/datasets/pancreas.py",
    "last_modified": "2022-12-04 12:06:43 -0500",
    "sha": "4bb8a7e04545a06c336d3d9364a1dd84fa2af1a4"
  },
  {
    "path": "tasks/denoising/datasets/pbmc.py",
    "last_modified": "2022-12-04 12:06:43 -0500",
    "sha": "4bb8a7e04545a06c336d3d9364a1dd84fa2af1a4"
  },
  {
    "path": "tasks/denoising/datasets/tabula_muris_senis.py",
    "last_modified": "2022-12-04 12:06:43 -0500",
    "sha": "4bb8a7e04545a06c336d3d9364a1dd84fa2af1a4"
  },
  {
    "path": "tasks/denoising/datasets/utils.py",
    "last_modified": "2022-11-15 17:19:16 -0500",
    "sha": "c2470ce02e6f196267cec1c554ba7ae389c0956a"
  },
  {
    "path": "tasks/denoising/methods/__init__.py",
    "last_modified": "2022-10-21 13:56:15 -0400",
    "sha": "b460ecb183328c857cbbf653488f522a4034a61c"
  },
  {
    "path": "tasks/denoising/methods/alra.R",
    "last_modified": "2022-05-16 15:10:42 -0400",
    "sha": "ba06cf71b564eb23823a662341055dc5ac2be231"
  },
  {
    "path": "tasks/denoising/methods/alra.py",
    "last_modified": "2022-07-25 12:29:34 -0400",
    "sha": "411a416150ecabce25e1f59bde422a029d0a8baa"
  },
  {
    "path": "tasks/denoising/methods/baseline.py",
    "last_modified": "2022-10-21 13:56:15 -0400",
    "sha": "b460ecb183328c857cbbf653488f522a4034a61c"
  },
  {
    "path": "tasks/denoising/methods/dca.py",
    "last_modified": "2022-12-01 15:38:21 -0500",
    "sha": "aa2253779e9aa9cd178f54ac0f3b6ba521ecd59f"
  },
  {
    "path": "tasks/denoising/methods/knn_smoothing.py",
    "last_modified": "2022-11-14 11:54:15 -0500",
    "sha": "bbecf4e9ad90007c2711394e7fbd8e49cbd3e4a1"
  },
  {
    "path": "tasks/denoising/methods/magic.py",
    "last_modified": "2022-11-14 11:57:35 -0500",
    "sha": "2af9a4918ed3370859f71774558068961f6d22c6"
  },
  {
    "path": "tasks/denoising/metrics/__init__.py",
    "last_modified": "2021-01-19 13:31:20 -0500",
    "sha": "8e0600c516c392fa747137415b6a93b8af0f61d8"
  },
  {
    "path": "tasks/denoising/metrics/mse.py",
    "last_modified": "2022-11-15 17:19:16 -0500",
    "sha": "c2470ce02e6f196267cec1c554ba7ae389c0956a"
  },
  {
    "path": "tasks/denoising/metrics/poisson.py",
    "last_modified": "2022-12-04 12:06:43 -0500",
    "sha": "4bb8a7e04545a06c336d3d9364a1dd84fa2af1a4"
  }
]
EOT
    

bin/viash run src/common/get_method_info/config.vsh.yaml -- \
    --input "src/denoising" \
     --output $OUTPUT_DIR/"method_info.json"