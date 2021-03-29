from . import utils

import anndata
import contextlib
import logging
import subprocess
import tempfile

log = logging.getLogger("openproblems")


def _write_h5ad_patch(self, *args, **kwargs):
    self.obs.columns = self.obs.columns.astype(str)
    self.var.columns = self.var.columns.astype(str)
    self._write_h5ad(*args, **kwargs)


@utils.temporary(version="0.5.0")
def patch_anndata():
    """Temporary fix for https://github.com/theislab/anndata/pull/457."""
    log.debug("Patching AnnData.write_h5ad")
    anndata.AnnData._write_h5ad = anndata.AnnData.write_h5ad
    anndata.AnnData.write_h5ad = _write_h5ad_patch


def _download_aftp(
    download_url,
    timeout=None,
    base_name="download",
    ext="tmp",
    use_wget_if_available=False,
):
    """Download FTP using wget.

    Anonymous ftp user circumvents Travis' FTP ban.
    """
    with tempfile.NamedTemporaryFile(
        suffix="." + ext, prefix=base_name, delete=False
    ) as tmp:
        tmp_path = tmp.name

    wget_command_list = [
        "wget",
        "--ftp-user=anonymous",
        "--ftp-password=anonymous@singlecellopenproblems.org",
        download_url,
        "-O",
        tmp_path,
        "--no-verbose",
        "--passive-ftp",
    ]
    if timeout:
        wget_command_list += ["-T", str(timeout)]

    log.debug("Running: %s" % (" ".join(wget_command_list)))
    subprocess.call(wget_command_list)
    return tmp_path


@contextlib.contextmanager
def patch_datacache():
    """Workaround for Travis CI's FTP ban."""
    log.debug("Patching datacache.download._download_to_temp_file")
    import datacache

    _download_to_temp_file = datacache.download._download_to_temp_file
    datacache.download._download_to_temp_file = _download_aftp
    try:
        yield None
    finally:
        datacache.download._download_to_temp_file = _download_to_temp_file
