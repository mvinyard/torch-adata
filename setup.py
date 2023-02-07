
__module_name__ = "setup.py"
__doc__ = """PYPI package distribution setup."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify package version: --------------------------------------------------
__version__ = "0.0.18"


# -- import packages: ----------------------------------------------------------
import setuptools
import re
import os
import sys


# -- run setup: ----------------------------------------------------------------
setuptools.setup(
    name="torch-adata",
    version="0.0.18",
    python_requires=">3.7.0",
    author="Michael E. Vinyard",
    author_email="mvinyard@g.harvard.edu",
    url=None,
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="torch-adata",
    packages=setuptools.find_packages(),
    install_requires=[
        "anndata>=0.8",
        "licorice_font>=0.0.3",
        "pytorch-lightning>=1.7.7",
        "torch>=1.12",        
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
