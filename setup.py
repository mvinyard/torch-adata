
__module_name__ = "setup.py"
__doc__ = """PYPI package distribution setup."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard.ai@gmail.com"])


# -- import packages: ----------------------------------------------------------
import setuptools


# -- fetch requirements packages: ---------------------------------------------
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('torch_adata/__version__.py') as v:
    exec(v.read())


# -- run setup: ---------------------------------------------------------------
setuptools.setup(
    name="torch_adata",
    version=__version__,
    python_requires=">3.9.0",
    author="Michael E. Vinyard",
    author_email="mvinyard.ai@gmail.com",
    url="https://github.com/mvinyard/torch-adata",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="torch-adata: bridging AnnData-based data models to torch.",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
