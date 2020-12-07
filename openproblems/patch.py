import anndata
import logging
from . import utils

log = logging.getLogger("openproblems")


def _write_h5ad_patch(self, *args, **kwargs):
    self.obs.columns = self.obs.columns.astype(str)
    self.var.columns = self.var.columns.astype(str)
    self._write_h5ad(*args, **kwargs)


@utils.temporary(version="0.3.1")
def patch_anndata():
    """Temporary fix for https://github.com/theislab/anndata/pull/457."""
    log.debug("Patching AnnData.write_h5ad")
    anndata.AnnData._write_h5ad = anndata.AnnData.write_h5ad
    anndata.AnnData.write_h5ad = _write_h5ad_patch
