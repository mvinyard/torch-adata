
__module_name__ = "_pad_time_resolved_dataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: -------------------------------------------------------
import numpy as np
import torch


# import local dependencies: ---------------------------------------------
from ._use_X import _use_X, _tensorize
from ._fetch_labels_from_obs import _fetch_labels_from_obs


# supporting functions: sampling -----------------------------------------
def _format_sampled(df, adata, use_key):
    
    if use_key in adata.obsm_keys():
        return torch.Tensor(adata[df.index].obsm[use_key])
    elif use_key in adata.obs_keys():
        return _fetch_labels_from_obs(adata[df.index], use_key, return_obj=False)
    elif use_key == "X":
        return _tensorize(adata.X)
    else:
        print("Must provide a valid key")
        
def _sample_dataset(adata, time_key, use_key):

    grouped = adata.obs.groupby(time_key)
    X_dict = grouped.apply(_format_sampled, adata, use_key).to_dict()
    return X_dict


# supporting functions: padding ------------------------------------------
def _pad_count(df, time_key="Time point"):
    
    count_df = df[time_key].value_counts().to_frame().reset_index()
    count_df.columns = [time_key, "count"]
    count_df["n_pad"] = count_df["count"].max() - count_df["count"]
    
    return count_df

def _randomly_pad_cells_at_timepoint(df, counted, time_key):

    t = df[time_key].unique()[0]
    t_count = counted.loc[counted[time_key] == t]["n_pad"].values[0]
    return np.random.choice(df.index, t_count)

def _pad_cells(adata, time_key):
    """randomly select cells that function as padding cells"""

    df = adata.obs.copy()
    pad_count_df = _pad_count(df, time_key)

    max_t = pad_count_df[pad_count_df["n_pad"] == 0][time_key].values[0]
    pad_dict = (
        df.groupby(time_key)
        .apply(_randomly_pad_cells_at_timepoint, pad_count_df, time_key)
        .to_dict()
    )

    del pad_dict[max_t]

    return pad_dict

def _apply_padding(data_dict, pad_cell_dict, adata, use_key):
    for key in data_dict.keys():
        if key in pad_cell_dict.keys():
            data_pad = _use_X(adata[pad_cell_dict[key]], use_key, return_obj=False)
            data_dict[key] = torch.vstack([data_dict[key], data_pad])
    return torch.stack(list(data_dict.values()))


# Main module function: --------------------------------------------------
def _pad_time_resolved_dataset(adata, time_key, use_key, obs_key):

    """samples dataset for padding"""
    
    pad_cell_dict = _pad_cells(adata, time_key)
    
    X_dict = _sample_dataset(adata, time_key, use_key)    
    X = _apply_padding(X_dict, pad_cell_dict, adata, use_key)
    
    if obs_key:
        y_dict = _sample_dataset(adata, time_key, obs_key)
        y = _apply_padding(y_dict, pad_cell_dict, adata, obs_key)
        
        return X, y
    else:
        return X