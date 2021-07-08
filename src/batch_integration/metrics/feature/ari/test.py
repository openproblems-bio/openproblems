from os import path
import subprocess
import pandas as pd

metric = 'ari.tsv'

print(">> Running script")
out = subprocess.check_output([
    "./ari",
    "--adata", 'mnn.h5ad',
    '--hvgs', '2000',
    "--output", metric
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric)

print(">> Check that score makes sense")
result = pd.read_table(metric)
assert result.shape == (1, 5)
ari = result.loc[0, 'value']
assert 0 < ari < 1
assert ari == 0.2459195865045752

print(">> All tests passed successfully")
