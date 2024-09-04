import subprocess

def test_with_duplicates(exec_path):
  # Create a temporary file with duplicate entries
  with open("test_with_duplicates.bib", mode="w") as file:
    file.write("@article{duplicate,\n author = {Duplicate, A.},\n title = {Duplicate article},\n year = {2022},\n}\n@article{duplicate,\n author = {Duplicate, A.},\n title = {Duplicate article},\n year = {2022},\n}\n")

  # Test basic functionality without merging
  result = subprocess.run(
    [exec_path, "--library", "test_with_duplicates.bib", "--library_v1", ""],
    capture_output=True,
    check=True,
    text=True
  )
  print(result.stdout, flush=True)
  
  assert "Read input bibtex file" in result.stdout, "Reading input failed"
  assert not "Merge with v1 library" in result.stdout, "Merging failed"
  assert "Remove duplicates" in result.stdout, "Duplicate removal failed"
  assert "Write to file" in result.stdout, "Writing output failed"

  # Check the output file to make sure duplicates are removed
  with open("test_with_duplicates.bib", "r") as f:
    contents = f.read()
    count = contents.count("article{")
    assert count == 1, f"Count should be 1 but is {count}"

def test_merge(exec_path, lib_v1_url):
  # Create a temporary file with duplicate entries
  with open("test_merge.bib", mode="w") as file:
    file.write("@article{entry,\n author = {Duplicate, A.},\n title = {Duplicate article},\n year = {2022},\n}\n")

  # Test basic functionality without merging
  result = subprocess.run(
    [exec_path, "--library", "test_merge.bib", "--library_v1", lib_v1_url],
    capture_output=True,
    check=True,
    text=True
  )
  print(result.stdout, flush=True)
  
  assert "Read input bibtex file" in result.stdout, "Reading input failed"
  assert "Merge with v1 library" in result.stdout, "Merging failed"
  assert "Remove duplicates" in result.stdout, "Duplicate removal failed"
  assert "Write to file" in result.stdout, "Writing output failed"

  # Check the output file to make sure duplicates are removed
  with open("test_merge.bib", "r") as f:
    contents = f.read()
    count = contents.count("article{")

    assert count > 1, f"Count should be greater than 1 but is {count}"

test_merge(
  exec_path=meta["executable"],
  lib_v1_url="https://raw.githubusercontent.com/openproblems-bio/openproblems/main/main.bib"
)
test_with_duplicates(
  exec_path=meta["executable"]
)

print("All tests passed!", flush=True)