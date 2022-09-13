from setuptools import setup
import re
import os
import sys


setup(
    name="torch-adata",
    version="0.0.12",
    python_requires=">3.6.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url=None,
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="torch-adata",
    packages=[
        "torch_adata",
        "torch_adata._functions",
    ],
    
    install_requires=[
        "anndata>=0.7.8",
        "torch>=1.12",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
