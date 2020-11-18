import anndata
import packaging.version
from .version import __version__


def _write_h5ad_patch(self, *args, **kwargs):
    self.obs.columns = self.obs.columns.astype(str)
    self.var.columns = self.var.columns.astype(str)
    self._write_h5ad(*args, **kwargs)


def patch_anndata():
    """Temporary fix for https://github.com/theislab/anndata/pull/457."""
    assert packaging.version.parse(__version__) < packaging.version.parse("0.3.1")
    anndata.AnnData._write_h5ad = anndata.AnnData.write_h5ad
    anndata.AnnData.write_h5ad = _write_h5ad_patch
