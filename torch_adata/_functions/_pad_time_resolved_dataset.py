
__module_name__ = "_pad_time_resolved_dataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import numpy as np
import torch

#### ------------------------------------------ ####

def _format_sampled(df, adata, use_key):
    if use_key in adata.obsm_keys():
        return torch.Tensor(adata[df.index].obsm[use_key])
    elif use_key in adata.obs.columns:
        return torch.Tensor(df[use_key])
    
def _sample_dataset(adata, time_key, data_key, weight_key): # , fate_bias_key):

    grouped = adata.obs.groupby(time_key)
    
    X_dict = grouped.apply(_format_sampled, adata, data_key).to_dict()
    W_dict = grouped.apply(_format_sampled, adata, weight_key).to_dict()
#     F_dict = grouped.apply(_format_sampled, adata, fate_bias_key).to_dict()

    return X_dict, W_dict # , F_dict

#### ------------------------------------------ ####

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

#### ------------------------------------------ ####


def _pad_time_resolved_dataset(adata, time_key, data_key, weight_key, fate_bias_key):

    """samples dataset for padding"""

    X_dict, W_dict = _sample_dataset(adata, time_key, data_key, weight_key) # , F_dict , fate_bias_key)
    pad_cell_dict = _pad_cells(adata, time_key)
        
    for key in X_dict.keys():
        if key in pad_cell_dict.keys():
            x_pad = torch.Tensor(adata[pad_cell_dict[key]].obsm[data_key])
            w_pad = torch.Tensor(adata[pad_cell_dict[key]].obs[weight_key])
#             f_pad = torch.Tensor(adata[pad_cell_dict[key]].uns[fate_bias_key])

            X_dict[key] = torch.vstack([X_dict[key], x_pad])
            W_dict[key] = torch.hstack([W_dict[key], w_pad])
#             if float(key) == 2.0:
#                 F_dict[key] = torch.hstack([F_dict[key], f_pad])
                
            
    X = torch.stack(list(X_dict.values()))
    W = torch.stack(list(W_dict.values()))
#     F = torch.Tensor(list(F_dict.values())[0])

    return X, W # , F