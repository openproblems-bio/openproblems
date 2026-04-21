from pathlib import Path
import subprocess

from openproblems.project import read_task_metadata, render_task_readme_qmd

## VIASH START
par = {
    "input": "path/to/input",
    "output": "path/to/input/README.md",
}
## VIASH END

print("Read task metadata", flush=True)
metadata = read_task_metadata(par["input"])

print("Render README.qmd content", flush=True)
qmd_content = render_task_readme_qmd(metadata)

output_path = Path(par["output"])
output_path.parent.mkdir(parents=True, exist_ok=True)
qmd_path = output_path.with_suffix(".qmd")

print("Write README.qmd to file", flush=True)
qmd_path.write_text(qmd_content)

print("Render README.qmd to README.md", flush=True)
out = subprocess.run(
    ["quarto", "render", str(qmd_path), "--output", "-"],
    capture_output=True,
    text=True,
)

if out.returncode != 0:
    stdout = (out.stdout or "").strip()
    stderr = (out.stderr or "").strip()
    details = []
    if stdout:
        details.append(f"stdout:\n{stdout}")
    if stderr:
        details.append(f"stderr:\n{stderr}")
    if not details:
        details.append("No output captured from Quarto")
    raise RuntimeError(
        "Failed to render README.qmd with Quarto.\n\n" + "\n\n".join(details)
    )

output_path.write_text(out.stdout)
