from os import path
import subprocess
import pandas as pd
import numpy as np
import scanpy as sc

np.random.seed(42)

metric = 'cell_cycle_conservation'
metric_file = metric + '.tsv'

print(sc.read('combat.h5ad'))

print(">> Running script")
out = subprocess.check_output([
    "./" + metric,
    "--adata", 'combat.h5ad',
    "--organism", "human",
    "--output", metric_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric_file)

print(">> Check that score makes sense")
result = pd.read_table(metric_file)
assert result.shape == (1, 4)
score = result.loc[0, 'value']
print(score)

assert 0 <= score <= 1


print(">> All tests passed successfully")
