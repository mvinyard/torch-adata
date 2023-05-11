
__module_name__ = "setup.py"
__doc__ = """PYPI package distribution setup."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: ----------------------------------------------------------
import setuptools
import re
import os
import sys


# -- run setup: ----------------------------------------------------------------
setuptools.setup(
    name="torch-adata",
    version="0.0.23rc0",
    python_requires=">3.9.0",
    author="Michael E. Vinyard",
    author_email="mvinyard@g.harvard.edu",
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
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
