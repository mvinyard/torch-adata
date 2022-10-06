# Ancilliary functions


```python
import torch_adata

idx = torch_adata.idx(dataset)
```


```python
batch = torch_adata.dummy_batch(dataset)
```

```python
# example
train, val = torch_adata.split(dataset)
```
