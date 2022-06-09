from setuptools import find_packages
from setuptools import setup

import os

install_requires = [
    "numpy>=1.19,<1.23",
    "scikit-learn>=0.24,<1.2",
    "anndata>=0.7.5,<0.8",
    "scprep>=1.2.0",
    "scipy>=1.5,<1.9",
    "scanpy>=1.6",
    "louvain==0.7.*",
    "decorator<5.0",
    "memory-profiler==0.60",
    "colorama==0.4.*",
    "packaging==21.3",
    "umap-learn==0.5.*",
]

r_requires = [
    "rpy2<3.4.3",
    "anndata2ri==1.0.6",
]

evaluate_requires = ["snakemake==7.8.*"]

process_requires = ["numpyencoder==0.3.*"]

test_requires = [
    "pytest==7.1.*",
    "pytest-cov==3.0.*",
    "black==22.3.0",
    "coverage==6.4.*",
    "codecov==2.1.*",
    "parameterized==0.8.*",
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
    python_requires=">=3.6",
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
