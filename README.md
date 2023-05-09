# ![torch-adata-logo](/docs/imgs/torch-adata.logo.github.png)

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/torch-adata.svg)](https://pypi.python.org/pypi/torch-adata/)
[![PyPI version](https://badge.fury.io/py/torch-adata.svg)](https://badge.fury.io/py/torch-adata)
[![Documentation Status](https://readthedocs.org/projects/torch-adata/badge/?version=latest)](https://torch-adata.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Create [`PyTorch Datasets`](https://pytorch.org/tutorials/beginner/basics/data_tutorial.html) from [`AnnData`](https://anndata.readthedocs.io/en/latest/)

## Installation

Install from PYPI (current version: **[`0.0.23`](https://pypi.org/project/torch-adata/)**):
```BASH
pip install torch-adata
```

Install the developer version:
```BASH
git clone https://github.com/mvinyard/torch-adata.git; cd torch-adata;
pip install -e .
```

## The main API

The primary class is the [`AnnDataset`](https://github.com/mvinyard/torch-adata/blob/main/torch_adata/_core/_AnnDataset.py). This is a subclass of the widely-used [`torch.utils.data.Dataset`](https://pytorch.org/tutorials/beginner/basics/data_tutorial.html). The PyTorch `Dataset` module enables us to take advantage of built-in multiprocessing and other organizational tricks that ultimately standardize workflows and enable reproducibility.

![torch-adata-concept-overview](/docs/imgs/torch-adata.concept_overview.png)

```python
import anndata as a
import torch_adata

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = torch_adata.AnnDataset(adata, use_key="X_pca", groupby="time", obs_keys=["affinity"])
```
```
[ torch-adata ]: AnnDataset object with 7131 samples
----------------------------------------------------
Grouped by: 'time' with attributes:
 - X (use_key = 'X_pca') torch.Size([3, 7131, 50])
 - obs: affinity: torch.Size([3, 7131, 1])
```

#
There is an additional approach to this dubbed [`AnnLoader`](https://github.com/scverse/anndata/blob/master/anndata/experimental/pytorch/_annloader.py), highlighted by [Sergei Rybakov](https://github.com/koncopd) in [Interfacing pytorch models with anndata](https://anndata-tutorials.readthedocs.io/en/latest/annloader.html)


**For more information, please visit the [documentation](https://torch-adata.readthedocs.io/en/latest/index.html)!**

**Problem?** Open an [issue](https://github.com/mvinyard/torch-adata/issues/new)
