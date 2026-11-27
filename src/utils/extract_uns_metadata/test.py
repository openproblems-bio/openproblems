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

def test_run_zarr(run_component, file_raw, tmp_path):
  import anndata as ad

  zarr_path = tmp_path / "dataset.anndata.zarr"
  ad.read_h5ad(input_path).write_zarr(zarr_path)

  schema = yaml.safe_load(file_raw.read_text())
  schema["info"]["format"]["type"] = "anndata_zarr"
  file_raw.write_text(yaml.dump(schema))

  output_path = tmp_path / "meta.yaml"

  run_component([
    "--input", str(zarr_path),
    "--schema", str(file_raw),
    "--output", str(output_path),
  ])

  assert output_path.exists(), "Output path does not exist"

  with open(output_path, "r") as f:
    out = yaml.safe_load(f)
  assert out["uns"]["file_size"] > 0, "File size of a zarr store should be the sum of its files"

if __name__ == "__main__":
  sys.exit(pytest.main([__file__]))
