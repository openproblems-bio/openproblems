import anndata
import json
import sys

import openproblems
import openproblems.test.utils


def evaluate_metric(adata, metric):
    result = float(metric(adata))
    return result


def main():
    openproblems.test.utils.ignore_numba_warnings()

    metric_name = sys.argv[0]
    adata_file = sys.argv[1]
    output_file = sys.argv[2]

    metric = eval(metric_name)
    adata = anndata.read_h5ad(adata_file)

    result = evaluate_metric(adata, metric)

    with open(output_file, "w") as handle:
        json.dump(result, handle, indent=4)


if __name__ == "__main__":
    main()
