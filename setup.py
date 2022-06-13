from setuptools import find_packages
from setuptools import setup

import os

install_requires = [
    "numpy>=1.17.0",
    "scikit-learn>=0.19.1",
    "anndata>=0.7.6,<0.8",
    "scprep>=1.2.0",
    "scipy",
    "scanpy>=1.6",
    "louvain>=0.7",
    "decorator<5.0",
    "memory-profiler",
    "colorama>=0.3.9",
    "packaging",
    "umap-learn>=0.5.1",
]

r_requires = [
    "rpy2<3.4.3",
    "anndata2ri>=1.0.6",
]

evaluate_requires = ["snakemake"]

process_requires = ["numpyencoder"]

test_requires = [
    "pytest",
    "pytest-cov",
    "black==22.3.0",
    "coverage",
    "codecov",
    "parameterized>=0.7.4",
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
        "evaluate": evaluate_requires,
        "process": process_requires,
    },
    entry_points={
        "console_scripts": ["openproblems-cli=openproblems.api.main:main"],
    },
    test_suite="nose2.collector.collector",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords=[
        "computational-biology",
    ],
)
