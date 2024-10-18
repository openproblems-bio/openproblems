import sys
import pytest
import yaml

## VIASH START
## VIASH END

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"

@pytest.fixture
def file_raw(tmp_path):
  file_raw_content = {
    "type": "file",
    "label": "Raw dataset",
    "summary": "An unprocessed dataset as output by a dataset loader.",
    "description": "This dataset contains raw counts and metadata as output by a dataset loader.",
    "info": {
      "format": {
        "type": "h5ad",
        "layers": [
          {
            "type": "integer",
            "name": "counts",
            "description": "Raw counts",
            "required": True
          }
        ],
        "obs": [
          {
            "type": "string",
            "name": "celltype",
            "description": "Classification of the cell type based on its characteristics and function within the tissue or organism.",
            "required": True
          }
        ]
      }
    }
  }
  file_raw_path = tmp_path / "file_raw.yaml"
  with open(file_raw_path, "w") as f:
    f.write(yaml.dump(file_raw_content))

  return file_raw_path

def test_run(run_component, file_raw, tmp_path):
  output_path = tmp_path / "meta.yaml"

  run_component([
    "--input", input_path,
    "--schema", str(file_raw),
    "--output", str(output_path),
  ])

  assert output_path.exists(), "Output path does not exist"

if __name__ == "__main__":
  sys.exit(pytest.main([__file__]))
