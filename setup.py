
__module_name__ = "setup.py"
__doc__ = """PYPI package distribution setup."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard.ai@gmail.com"])


# -- import packages: ----------------------------------------------------------
import setuptools
import re
import os
import sys


# -- run setup: ----------------------------------------------------------------
setuptools.setup(
    name="torch-adata",
    version="0.0.24rc1",
    python_requires=">3.9.0",
    author="Michael E. Vinyard",
    author_email="mvinyard.ai@gmail.com",
    url=None,
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="torch-adata",
    packages=setuptools.find_packages(),
    install_requires=[
        "anndata>=0.9.1",
        "licorice_font>=0.0.3",
        "lightning>=2.0.1",
        "torch>=2.0",
#         "numpy==1.23",
        "scanpy==1.9.3",
        "scikit-learn==1.2.2",
        "webfiles",
	"vinplots",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
