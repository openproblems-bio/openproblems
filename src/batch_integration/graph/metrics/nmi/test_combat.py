from os import path
import subprocess
import pandas as pd
import numpy as np

np.random.seed(42)

metric = 'nmi'
metric_file = metric + '.tsv'

print(">> Running script")
out = subprocess.check_output([
    "./" + metric,
    "--input", 'combat.h5ad',
    "--output", metric_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric_file)
result = pd.read_table(metric_file)

print(">> Check that score makes sense")
assert result.shape == (1, 4)
score = result.loc[0, 'value']
print(score)

assert 0 <= score <= 1

print(">> All tests passed successfully")
