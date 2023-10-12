import sys
import re
import pytest
import json
import subprocess

## VIASH START
## VIASH END

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"

@pytest.fixture
def schema(tmp_path):
  schema = tmp_path / "schema.yaml"
  schema.write_text("""
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  label: "Preprocessed dataset"
  slots:
    layers: 
      - type: integer
        name: counts
        description: Raw counts
        required: true
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
        required: true
""")
  return schema

@pytest.fixture
def error_schema(tmp_path):
  schema = tmp_path / "schema.yaml"
  schema.write_text("""
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  label: "Preprocessed dataset"
  slots:
    X:
      type: double
      description: Normalized expression values
      required: true
    layers: 
      - type: integer
        name: counts
        description: Raw counts
        required: true
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
        required: true
      - type: string
        name: error_test
        description: "A made up uns variable to test if error is picked up"
        required: true
  """)
  return schema

def test_run(run_component, tmp_path, schema):
  output_path = tmp_path / "output.h5ad"

  run_component([
    "--input", input_path,
    "--schema", str(schema),
    "--output", str(output_path)
  ])

  assert output_path.exists(), "Output path does not exist"

def test_error(run_component, tmp_path, error_schema):
  output_checks = tmp_path / "checks.json"
  output_path = tmp_path / "output.h5ad"

  with pytest.raises(subprocess.CalledProcessError) as err:
    run_component([
      "--input", input_path,
      "--schema", str(error_schema),
      "--stop_on_error", "true",
      "--checks", str(output_checks),
      "--output", str(output_path)
    ])
    assert err.value.exitcode > 0

  assert output_checks.exists(), "Output checks file does not exist"
  assert not output_path.exists(), "Output path does not exist"

  with open(output_checks, "r") as f:
      out = json.load(f)
      assert out["exit_code"] > 0
      assert out["data_schema"] == "not ok"


if __name__ == "__main__":
  sys.exit(pytest.main([__file__]))
