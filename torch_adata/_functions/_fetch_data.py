
__module_name__ = "_fetch_data.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages --------------------------------------------------------
import torch
import numpy
import pandas
from sklearn.preprocessing import OneHotEncoder


def _one_hot_encode_labels(y):
    """first format for one-hot encoding then actually do the OHE"""

    OneHot = OneHotEncoder()
    y_one_hot = OneHot.fit_transform(y).toarray()
    return OneHot, torch.Tensor(y_one_hot)


# fetch labels from obs (y) ----------------------------------------------
def _fetch_y_from_obs(adata, obs_key, onehot_encode_y=False, return_obj=True):
    
    """
    Fetch data labels (y) from adata.obs[obs_key]
    
    Parameters:
    -----------
    adata
        AnnData object
        type: anndata._core.anndata.AnnData
    
    obs_key
        pandas DataFrame column accessor in adata.obs_keys()
        type: str
    
    Returns:
    --------
    y
        torch.Tensor
    
    Notes:
    ------
    (1) Creating this function as a wrapper of OHE because I'm not sure if
        there may be other uses beyond OHE
    """
    
    y = _to_np_array(adata.obs[obs_key].values).reshape(-1, 1)
    
    if onehot_encode_y:
        OneHot, y_one_hot = _one_hot_encode_labels(y)
        if return_obj:
            return OneHot, y_one_hot
    else:
        if return_obj:
            return None, torch.Tensor(y)
        else:
            return torch.Tensor(y)
        

# fetch data (X) ---------------------------------------------------------
def _is_numpy_array(x):
    return x.__class__ is numpy.ndarray

def _is_pandas_categorical(x):
    x.__class__ is pandas.Categorical

def _to_np_array(x):
    if _is_numpy_array(x):
        return x
    elif _is_pandas_categorical(x):
        return x.to_numpy()
    else:
        return x.toarray()


def _tensorize(x):
    return torch.Tensor(_to_np_array(x))


def _use_X(adata, use_key="X", onehot_encode_y=False, return_obj=False):
    
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
        return _fetch_y_from_obs(adata, use_key, onehot_encode_y=False, return_obj=False)

    else:
        print("No suitable array found for dataset...")
