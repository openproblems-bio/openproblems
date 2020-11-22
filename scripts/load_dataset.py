import openproblems
import sys


def main(task_name, dataset_name, output_file):
    """Load and save a dataset to disk."""
    openproblems.data.no_cleanup()
    task = eval("openproblems.tasks.{}".format(task_name))
    datasets = getattr(task, "datasets")
    dataset = getattr(datasets, dataset_name)
    adata = dataset(test=False)
    adata.write_h5ad(output_file)


if __name__ == "__main__":
    main(*sys.argv[1:])
