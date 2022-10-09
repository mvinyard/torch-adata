
__module_name__ = "_fetch_data.py"
__doc__ = """Data fetching module. Key module called within torch_adata.AnnDataset.__init__()."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import torch
import pandas
import anndata


# -- import local dependencies: ----------------------------------------------------------
from ._group_and_pad_adata import group_and_pad_adata
from ._one_hot_encode import one_hot_encode
from ._data_typing import tensorize, to_np_array


# -- fetch X from adata: -----------------------------------------------------------------
def fetch_X(adata, use_key):
    """
    Flexibly fetch a data matrix to use as "X" from adata.

    Parameters:
    -----------
    adata

    use_key
        type: str

    Returns:
    --------
    X
        type: torch.Tensor

    Notes:
    ------
    (1) Suggested: "X_pca" (or similar).
    (2) For large datasets, using "X" can be very slow with a high memory requirement.
    """

    if use_key == "X":
        return tensorize(adata.X)
    if use_key in adata.layers:
        return tensorize(adata.layers[use_key])
    if use_key in adata.obsm_keys():
        return tensorize(adata.obsm[use_key])
    if use_key in adata.obs_keys():
        return fetch_from_obs(adata, use_key, one_hot=False)
    print("No suitable array found for dataset...")


# -- fetch from adata.obs table: ----------------------------------------------------------
def fetch_from_obs(
    obs_df: pandas.DataFrame, obs_key: str, one_hot: bool = False
) -> torch.Tensor:

    """
    Fetch a column from adata.obs.

    Parameters:
    -----------
    obs_df
        type: pandas.DataFrame

    obs_key
        type: str

    one_hot
        type: bool
        default: False

    Returns:
    --------
    y
        type: torch.Tensor
    """

    y = to_np_array(obs_df[obs_key].values).reshape(-1, 1)

    if y.dtype == "O" or one_hot:
        return one_hot_encode(y)[1]
    return torch.Tensor(y)


def fetch_from_multiple_obs_keys(adata, obs_keys, attr_names, one_hot):

    """
    Fetch obs data from multiple passed obs_keys.

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    obs_keys
        type: list([str, ..., str])

    attr_names:
        type: list([str, ..., str])

    one_hot_encode:
        type: list([bool, ..., bool])

    Returns:
    --------
    X, obs_stacked
    """

    obs_data = {}
    for n, key in enumerate(obs_keys):
        val = fetch_from_obs(adata.obs, obs_key=key, one_hot=one_hot[n])
        obs_data[attr_names[n]] = val

    return obs_data

def _subset_obs_names(attr_names, obs_keys)->list([str, ..., str]):
    return attr_names[:len(obs_keys)]

def _isolate_aux_names(attr_names, obs_keys) -> list([str, ..., str]):
    if len(attr_names) > len(obs_keys):
        return [attr_names[i] for i in range(len(obs_keys), len(attr_names))]

def fetch_from_grouped_adata(
    adata: anndata.AnnData,
    groupby: str,
    use_key: str,
    obs_keys: list([str, ..., str]),
    aux_keys: list([str, ..., str]),
    attr_names: list([str, ..., str]),
    one_hot: list([bool, ..., bool]),
):
    """
    Fetch both X and ancilliary obs data from a grouped anndata object.

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    groupby
        type: str

    use_key
        type: str

    obs_keys:
        type: list([str, ..., str])

    attr_names:
        obs_keys listed first, then aux_keys.
        type: list([str, ..., str])

    one_hot_encode:
        type: list([bool, ..., bool])

    Returns:
    --------
    X, obs_stacked
    """
    # -- (1) group adata.obs by key, then pad dataset indices: ---------------------------
    idx_dict = group_and_pad_adata(adata, groupby)
    
    # -- (2) create receptacles for stacked data: ----------------------------------------
    
    print("attr names", attr_names)
    obs_names = _subset_obs_names(attr_names, obs_keys)
    aux_names = _isolate_aux_names(attr_names, obs_keys)
    print("obs names", obs_names)
    print("aux names", aux_names)
    
    X_list = []
    if obs_names:
        obs_dict = {key: [] for key in obs_names}
#     aux_list = []
    if aux_names:
        aux_dict = {key: [] for key in aux_names}

    # -- (3) grab X and each obs vector for each group: ----------------------------------
    for group_idx in idx_dict.values():
        group_adata = adata[group_idx]
        X_list.append(fetch_X(group_adata, use_key))
        if obs_keys:
            obs_data = fetch_from_multiple_obs_keys(
                group_adata, obs_keys, attr_names, one_hot
            )
            for key in obs_dict.keys():
                obs_dict[key].append(obs_data[key])
        if aux_keys:
            for aux_attr in aux_names:
                aux_dict[aux_attr].append(fetch_X(group_adata, aux_attr))

    # -- (4) restack the data you just grabbed: ------------------------------------------
    X = torch.stack(X_list)
    if obs_keys and aux_keys:
        obs_stacked = {key: torch.stack(obs_dict[key]) for key in obs_names}
        aux_stacked = {key: torch.stack(aux_dict[key]) for key in aux_names}
        return X, obs_stacked, aux_stacked
    if (obs_keys) and (not aux_keys):
        obs_stacked = {key: torch.stack(obs_dict[key]) for key in obs_names}
        return X, obs_stacked, None
    if (aux_keys) and (not obs_keys):
        aux_stacked = {key: torch.stack(aux_dict[key]) for key in aux_names}
        return X, None, aux_stacked
    else:
        return X, None, None


# -- fetch controller class: -------------------------------------------------------------
class Fetch:
    """Controller class for various fetching functions."""
    def __init__(self, adata):

        """
        Parameters:
        -----------
        adata
            type: anndata.AnnData

        Returns:
        --------
        None
        """

        self.adata = adata
        self.obs_df = adata.obs

    def X(self, use_key):
        """Fetch X"""
        return fetch_X(self.adata, use_key)

    def obs(self, obs_key: str, one_hot: bool = False):
        """Fetch a col from obs"""
        return fetch_from_obs(self.obs_df, obs_key, one_hot)

    def multi_obs(self, obs_keys, attr_names, one_hot):
        """Fetch multiple cols from obs"""
        return fetch_from_multiple_obs_keys(
            self.adata, obs_keys, attr_names, one_hot
        )

    def grouped_adata(self, groupby, use_key, obs_keys, aux_keys, attr_names, one_hot):
        """
        Fetch X and/or multiple cols from obs using a groupby
        accessor for adata.obs
        """
        return fetch_from_grouped_adata(
            self.adata, groupby, use_key, obs_keys, aux_keys, attr_names, one_hot
        )

    def update_obs_attrs(self, dataset, obs_data):
        """
        Add attributes to the dataset from obs_keys to be returned
        upon calling dataset.__getitem__()
        """
        for attr_name, val in obs_data.items():
            setattr(dataset, attr_name, val)
    
    def update_aux_attrs(self, dataset, aux_data):
        """
        Add attributes to the dataset from aux_keys to be returned
        upon calling dataset.__getitem__()
        """
        for attr_name, val in aux_data.items():
            setattr(dataset, attr_name, val)
            
    def update_attrs(self, dataset, obs_data, aux_data):
        if obs_data:
            self.update_obs_attrs(dataset, obs_data)
        if aux_data:
            self.update_aux_attrs(dataset, aux_data)
            