import os
import sys
from setuptools import setup, find_packages

install_requires = [
    "numpy>=1.17.0",
    "scikit-learn>=0.19.1",
    "anndata",
    "mnnpy",
    "scprep",
    "scipy",
    "scanpy",
    "decorator",
    "rpy2",
    "scIB @ git+https://github.com/theislab/scib@master",
    "harmonicalignment @ git+https://github.com/KrishnaswamyLab/harmonic-alignment#subdirectory=python",
]

evaluate_requires = ["tabulate"]

test_requires = ["nose2", "parameterized", "black"]

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
        "test": test_requires,
        "doc": doc_requires,
        "evaluate": evaluate_requires,
    },
    test_suite="nose2.collector.collector",
    long_description=readme,
    long_description_content_type="text/markdown",
    keywords=["computational-biology",],
)
