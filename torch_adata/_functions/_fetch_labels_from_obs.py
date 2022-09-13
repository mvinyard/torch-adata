
__module_name__ = "_fetch_labels_from_obs.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages --------------------------------------------------------
from sklearn.preprocessing import OneHotEncoder
import torch


def _one_hot_encode_labels(adata, obs_key):

    """first format for one-hot encoding then actually do the OHE"""

    X = adata.obs[obs_key].values.to_numpy().reshape(-1, 1)
    OneHot = OneHotEncoder()
    X_one_hot = OneHot.fit_transform(X).toarray()
    return OneHot, torch.Tensor(X_one_hot)


# fetch labels from obs (y) ----------------------------------------------
def _fetch_labels_from_obs(adata, obs_key, return_obj=True):
    
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
    
    
    
    OneHot, X_one_hot = _one_hot_encode_labels(adata, obs_key)
    
    if return_obj:
        return OneHot, X_one_hot
    else:
        return X_one_hot