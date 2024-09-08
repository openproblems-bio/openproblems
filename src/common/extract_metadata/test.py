import sys
import re
import pytest
import json
import subprocess

## VIASH START
## VIASH END

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"
schema_path = meta["resources_dir"] + "/file_raw.yaml"

def test_run(run_component, tmp_path):
  output_path = tmp_path / "meta.yaml"

  run_component([
    "--input", input_path,
    "--schema", schema_path,
    "--output", str(output_path),
  ])

  assert output_path.exists(), "Output path does not exist"


if __name__ == "__main__":
  sys.exit(pytest.main([__file__]))
