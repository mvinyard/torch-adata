# ![torch-adata-logo](/docs/imgs/torch-adata.logo.large.svg)

### *Create pytorch Datasets from* [`AnnData`](https://anndata.readthedocs.io/en/latest/)

<a href="https://github.com/mvinyard/torch-adata/"><img src="/docs/imgs/torch-adata.concept_overview.svg" alt="torch-adata-concept-overview" height="250" /></a>
## Example use of the base class

The base class, `AnnDataset` is a subclass of the widely-used `torch.utils.data.Dataset`. 

```python
import anndata as a
import torch_adata

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = torch_adata.AnnDataset(adata)
```

Returns sampled data `X_batch` as a `torch.Tensor`.
```python
# create a dummy index
idx = np.random.choice(range(dataset.__len__()), 5)
X_batch = dataset.__getitem__(idx)
```

#### `TimeResolvedAnnDataset`

Specialized class for time-resolved datasets. A subclass of the class, `AnnDataset`.

```python
import anndata as a
import torch_adata as ta

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = torch_adata.TimeResolvedAnnDataset(adata, time_key="Time point")
```

## Installation

Install from PYPI:
```BASH
pip install torch-adata
```

Install the developer version:
```BASH
git clone https://github.com/mvinyard/torch-adata.git; cd torch-adata;
pip install -e .
```
