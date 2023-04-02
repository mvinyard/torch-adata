
__module_name__ = "_fetch_data.py"
__doc__ = """Data fetching module. Key module called within torch_adata.AnnDataset.__init__()."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import torch
import anndata


# -- import local dependencies: ----------------------------------------------------------
from .._core._core_ancilliary._fetch_data import Fetch


# -- fetch X from adata: -----------------------------------------------------------------
def fetch(adata: anndata.AnnData, use_key: str)->torch.Tensor:
    
    f = Fetch(adata)
    return f.X(use_key)
            