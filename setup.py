import os
from setuptools import setup, find_packages

install_requires = [
    "numpy>=1.17.0",
    "scikit-learn>=0.19.1",
    "anndata",
    "scprep>=1.0.10",
    "scipy",
    "scanpy>=1.6",
    "decorator",
    "memory-profiler",
    "parameterized",
]

r_requires = [
    "rpy2",
    "scIB @ git+https://github.com/theislab/scib@master",
    "anndata2ri>=1.0.4",
]

evaluate_requires = ["snakemake"]

test_requires = [
    "nose2",
    "black",
    "coverage",
    "coveralls",
]

doc_requires = [
    "sphinx>=2.2,<2.4",
    "sphinxcontrib-napoleon",
    "autodocsumm",
    "ipykernel",
    "nbsphinx",
]

version_py = os.path.join(os.path.dirname(__file__), "openproblems", "version.py")
version = open(version_py).read().strip().split("=")[-1].replace('"', "").strip()

readme = open("README.md").read()

setup(
    name="openproblems",
    version=version,
    packages=find_packages(),
    license="MIT",
    install_requires=install_requires,
    python_requires=">=3.5",
    extras_require={
        "test": test_requires + r_requires,
        "doc": doc_requires,
        "r": r_requires,
        "evaluate": evaluate_requires + r_requires,
    },
    test_suite="nose2.collector.collector",
    long_description=readme,
    long_description_content_type="text/markdown",
    keywords=[
        "computational-biology",
    ],
)
