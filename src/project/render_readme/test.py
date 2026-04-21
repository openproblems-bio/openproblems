import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

## VIASH START
meta = {
    "executable": "foo",
}
## VIASH END

task_template_repo = "https://github.com/openproblems-bio/task_template.git"

with TemporaryDirectory() as tmpdir:
    tmpdir_path = Path(tmpdir)
    task_template_path = tmpdir_path / "task_template"
    input_path = task_template_path / "src" / "api"
    output_path = Path(tmpdir) / "README.md"
    qmd_path = Path(tmpdir) / "README.qmd"

    clone_cmd = [
        "git",
        "clone",
        "--depth",
        "1",
        task_template_repo,
        str(task_template_path),
    ]

    print(">> Cloning task_template", flush=True)
    clone_out = subprocess.run(clone_cmd, capture_output=True, text=True)

    if clone_out.stdout:
        print(clone_out.stdout)
    if clone_out.stderr:
        print(clone_out.stderr)

    if clone_out.returncode:
        print(f"script: '{' '.join(clone_cmd)}' exited with an error.")
        exit(clone_out.returncode)

    cmd = [
        meta["executable"],
        "--input",
        str(input_path),
        "--output",
        str(output_path),
    ]

    print(">> Running the script as test", flush=True)
    out = subprocess.run(cmd, capture_output=True, text=True)

    if out.stdout:
        print(out.stdout)
    if out.stderr:
        print(out.stderr)

    if out.returncode:
        print(f"script: '{' '.join(cmd)}' exited with an error.")
        exit(out.returncode)

    print(">> Checking whether output files exist", flush=True)
    assert output_path.exists(), "README.md was not generated"
    assert qmd_path.exists(), "README.qmd was not generated"

    print(">> Checking file contents", flush=True)
    output_lines = output_path.read_text().splitlines()
    qmd_lines = qmd_path.read_text().splitlines()

    assert any("## Description" in line for line in output_lines), "README.md is missing description section"
    assert any("flowchart TB" in line for line in output_lines), "README.md is missing task graph"
    assert any("## File format:" in line for line in output_lines), "README.md is missing file format section"
    assert any("format: gfm" in line for line in qmd_lines), "README.qmd header is missing gfm format"

print("All checks succeeded!", flush=True)
