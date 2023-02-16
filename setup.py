from setuptools import find_packages
from setuptools import setup

import os

install_requires = [
    "anndata==0.8.*",
    "colorama==0.4.*",
    "decorator<5.0",  # pinned in #324
    "louvain==0.8.*",
    "matplotlib<3.7.0",
    "memory-profiler==0.60",
    "numpy>=1.21,<1.24",
    "packaging==21.3",
    "pandas==1.3.5",
    "python-igraph==0.10.*",
    "scanpy>=1.6",
    "scipy>=1.7,<1.10",
    "scikit-learn>=1.0,<1.2",
    "scprep>=1.2.2",
    "umap-learn==0.5.*",
    "requests==2.28.*",
]

r_requires = [
    "rpy2>=3.4,<3.4.3",
    "anndata2ri==1.0.6",
]

evaluate_requires = ["snakemake>=7.8,<7.17", "tabulate<0.9"]

process_requires = ["numpyencoder==0.3.*"]

test_requires = [
    "pytest==7.1.*",
    "pytest-cov>=3.0,<4.1",
    "pytest-timestamper==0.0.9",
    "black==22.10.0",
    "coverage>=6.4,<6.6",
    "codecov==2.1.*",
    "parameterized==0.8.*",
    "requests==2.28.*",
    "bibtexparser==1.4.*",
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
    python_requires=">=3.7",
    extras_require={
        "test": test_requires + r_requires,
        "r": r_requires,
        "evaluate": evaluate_requires,
        "process": process_requires,
    },
    entry_points={
        "console_scripts": ["openproblems-cli=openproblems.api.main:main"],
    },
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords=[
        "computational-biology",
    ],
)
