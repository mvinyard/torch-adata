# torch-adata
Create pytorch Datasets from AnnData

## Example use of the base class

The base class, `AnnDataset` is a subclass of the widely-used `torch.utils.data.Dataset`. The outputs of all `AnnDataset` classes and subclasses are designed to be directly compatible with the `torch.utils.data.DataLoader` module.

```python
import anndata as a
import torch_adata as ta

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = ta.AnnDataset(adata)
```

Returns data (`X` as a `torch.Tensor`) and the `pandas.DataFrame`; `adata.obs`.
```python
# create a dummy index
idx = np.random.choice(range(dataset.__len__()), 5)
X, obs = dataset.__getitem__(idx)
```

## Specialized classes

#### `GroupedAnnDataset`

A subclass of the base class, `AnnDataset`.

```python
import anndata as a
import torch_adata as ta

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = ta.GroupedAnnDataset(adata, groupby="batch")
```

Returns data as a dictionary of data with values as `torch.Tensor` and keys as each `groupby` category and the sampled `adata.obs` is again returned as a  `pandas.DataFrame`.

```python
# create a dummy index
idx = np.random.choice(range(dataset.__len__()), 5)
X_dict, obs = dataset.__getitem__(idx)
```

#### `TimeResolvedAnnDataset`

A subclass of the class, `GroupedAnnDataset`.

```python
import anndata as a
import torch_adata as ta

adata = a.read_h5ad("/path/to/data.h5ad")
dataset = ta.TimeResolvedAnnDataset(adata, time_key="Time point")
```

Returns the initial datapoint, `X0` as a `torch.Tensor`, the entire sample of the dataset as a dictionary of data with values as `torch.Tensor` and keys as each timepoint indicated by the `time_key`. Sampled `adata.obs` is again returned as a  `pandas.DataFrame`.

```python
# create a dummy index
idx = np.random.choice(range(dataset.__len__()), 5)
X0, X_dict, t, obs = dataset.__getitem__(idx)
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
