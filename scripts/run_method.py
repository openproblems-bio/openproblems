import openproblems
import sys
import json
import anndata


def main(task_name, method_name, input_file, output_h5ad, output_json):
    """Apply a method to a single dataset."""
    openproblems.data.no_cleanup()
    task = eval("openproblems.tasks.{}".format(task_name))
    methods = getattr(task, "methods")
    method = getattr(methods, method_name)
    adata = anndata.read_h5ad(input_file)

    output = openproblems.tools.decorators.profile(method)(adata)
    adata = output["result"]
    adata.write_h5ad(output_h5ad)
    del adata

    result = dict()
    result["Name"] = method.metadata["method_name"]
    result["Paper"] = method.metadata["paper_name"]
    result["Paper URL"] = method.metadata["paper_url"]
    result["Year"] = method.metadata["paper_year"]
    result["Code"] = method.metadata["code_url"]
    result["Version"] = method.metadata["code_version"]
    result["Memory (GB)"] = float(output["memory_mb"] / 1024)
    result["Memory leaked (GB)"] = float(output["memory_leaked_mb"] / 1024)
    result["Runtime (min)"] = float(output["runtime_s"] / 60)

    with open(output_json, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main(*sys.argv[1:])
