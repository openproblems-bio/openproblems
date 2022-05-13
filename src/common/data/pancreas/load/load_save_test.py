import sys
sys.path.append("./")
import load_save
import fake_anndata
from os.path import exists
import unittest
import scanpy as sc
from unittest.mock import patch


URL = "https://fake.url"
OUTPUT = "/tmp/output.h5ad"


class LoadDataTest(unittest.TestCase):
    @patch('scprep.io.download.download_url')
    @patch('scanpy.read')
    def test_load_data(self, read_method, download_url_method):
        download_url_method.return_value = None
        read_method.return_value = fake_anndata.generate_fake_anndata()
        response = load_save.load_data(URL)
        self.assertEqual((3, 3), response.X.shape)


class SaveDataTest(unittest.TestCase):
    def test_save_data(self):
        load_save.save_data(fake_anndata.generate_fake_anndata(), "./output.h5ad")
        self.assertTrue(exists("./output.h5ad"), "File should exists")


if __name__ == "__main__":
    unittest.main()
