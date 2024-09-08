import warnings
from pathlib import Path
import anndata as ad
import h5py
from scipy.sparse import csr_matrix
from anndata.experimental import read_elem, sparse_dataset


def read_anndata(
    file: str,
    backed: bool = False,
    **kwargs
) -> ad.AnnData:
    """
    Read anndata file
    :param file: path to anndata file in h5ad format
    :param kwargs: AnnData parameter to group mapping
    """
    assert Path(file).exists(), f'File not found: {file}'
    
    f = h5py.File(file, 'r')
    kwargs = {x: x for x in f} if not kwargs else kwargs
    if len(f.keys()) == 0:
        return ad.AnnData()
    # check if keys are available
    for name, slot in kwargs.items():
        if slot not in f:
            warnings.warn(
                f'Cannot find "{slot}" for AnnData parameter `{name}` from "{file}"'
            )
    adata = read_partial(f, backed=backed, **kwargs)
    if not backed:
        f.close()
    
    return adata


def read_partial(
    group: h5py.Group,
    backed: bool = False,
    force_sparse_types: [str, list] = None,
    **kwargs
) -> ad.AnnData:
    """
    Partially read h5py groups
    :params group: file group
    :params force_sparse_types: encoding types to convert to sparse_dataset via csr_matrix
    :params backed: read sparse matrix as sparse_dataset
    :params **kwargs: dict of slot_name: slot, by default use all available slot for the h5py file
    :return: AnnData object
    """
    if force_sparse_types is None:
        force_sparse_types = []
    elif isinstance(force_sparse_types, str):
        force_sparse_types = [force_sparse_types]
    slots = {}
    if backed:
        print('Read as backed sparse matrix...')
    
    for slot_name, slot in kwargs.items():
        print(f'Read slot "{slot}", store as "{slot_name}"...')
        if slot not in group:
            warnings.warn(f'Slot "{slot}" not found, skip...')
            slots[slot_name] = None
        else:
            elem = group[slot]
            iospec = ad._io.specs.get_spec(elem)
            if iospec.encoding_type in ("csr_matrix", "csc_matrix") and backed:
                slots[slot_name] = sparse_dataset(elem)
            elif iospec.encoding_type in force_sparse_types:
                slots[slot_name] = csr_matrix(read_elem(elem))
                if backed:
                    slots[slot_name] = sparse_dataset(slots[slot_name])
            else:
                slots[slot_name] = read_elem(elem)
    return ad.AnnData(**slots)

