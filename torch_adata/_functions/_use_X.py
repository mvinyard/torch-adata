
__module_name__ = "_use_X.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages --------------------------------------------------------
import torch
import numpy

from ._fetch_labels_from_obs import _fetch_labels_from_obs

# fetch data (X) ---------------------------------------------------------
def _is_numpy_array(x):
    return x.__class__ is numpy.ndarray

def _toarray(x):
    if _is_numpy_array(x):
        return x
    else:
        return x.toarray()


def _tensorize(x):
    return torch.Tensor(_toarray(x))


def _use_X(adata, use_key="X", return_obj=False):
    
    """
    Return data from AnnData as Tensor.

    Parameters:
    -----------
    adata
        AnnData object
        type: anndata._core.anndata.AnnData

    use_key
        "X" or key in adata.obsm_keys(), adata.obs_keys(),
        type: str
        default: "X"

    Returns:
    --------
    X
        Formatted data matrix.
        type: torch.Tensor
    """

    if use_key == "X":
        return _tensorize(adata.X)

    elif use_key in adata.obsm_keys():
        return _tensorize(adata.obsm[use_key])
    
    elif use_key in adata.obs_keys():
        return _fetch_labels_from_obs(adata, use_key, return_obj)

    else:
        print("No suitable array found for dataset...")
