import anndata
import os


def _cache_name(tempdir, task, dataset, test=None, method=None):
    if not isinstance(task, str):
        task = task.__name__.split(".")[-1]
    if not isinstance(dataset, str):
        dataset = dataset.__name__
    if method is not None:
        if not isinstance(method, str):
            method = method.__name__
        return os.path.join(tempdir, "{}_{}_{}.h5ad".format(task, dataset, method))
    else:
        return os.path.join(tempdir, "{}_{}_{}.h5ad".format(task, dataset, test))


def load(tempdir, task, dataset, test=None, method=None, dependency="a related test"):
    """Load a cached h5ad file."""
    data_path = _cache_name(tempdir, task, dataset, test=test, method=method)
    assert os.path.isfile(data_path), "Intermediate file missing. Did {} fail?".format(
        dependency
    )
    return anndata.read_h5ad(data_path)


def save(adata, tempdir, task, dataset, test=None, method=None):
    """Save an AnnData object to h5ad."""
    data_path = _cache_name(tempdir, task, dataset, test=test, method=method)
    adata.write_h5ad(data_path)


def delete(tempdir, task, dataset, test=None, method=None):
    """Delete a cached AnnData object."""
    data_path = _cache_name(tempdir, task, dataset, test=test, method=method)
    try:
        os.remove(data_path)
    except FileNotFoundError:
        pass
