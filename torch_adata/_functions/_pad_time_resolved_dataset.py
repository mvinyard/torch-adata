
__module_name__ = "_pad_time_resolved_dataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import numpy as np
import torch


def _randomly_pad_cells_at_timepoint(df, counted, time_key):

    t = df[time_key].unique()[0]
    t_count = counted.loc[counted[time_key] == t]["n_pad"].values[0]
    return np.random.choice(df.index, t_count)


def _format_unsampled_data(df, adata, use_key):
    return torch.Tensor(adata[df.index].obsm[use_key])


def _sample_dataset(adata, time_key="Time point", use_key="X_pca"):

    X_dict = (
        adata.obs.groupby(time_key)
        .apply(_format_unsampled_data, adata, use_key)
        .to_dict()
    )

    return X_dict


def _pad_cells(adata, time_key):
    """randomly select cells that function as padding cells"""

    df = adata.obs.copy()

    count_df = df[time_key].value_counts().to_frame().reset_index()
    count_df.columns = [time_key, "count"]
    count_df["n_pad"] = count_df["count"].max() - count_df["count"]

    max_t = count_df[count_df["n_pad"] == 0][time_key].values[0]
    pad_dict = (
        df.groupby(time_key)
        .apply(_randomly_pad_cells_at_timepoint, count_df, time_key)
        .to_dict()
    )

    del pad_dict[max_t]

    return pad_dict


def _pad_time_resolved_dataset(adata, time_key="Time point", use_key="X_pca"):

    """samples dataset for padding"""

    X_dict = _sample_dataset(adata, time_key, use_key)
    pad_cell_dict = _pad_cells(adata, time_key)

    for key in X_dict.keys():
        if key in pad_cell_dict.keys():
            x_pad = torch.Tensor(adata[pad_cell_dict[key]].obsm[use_key])
            X_dict[key] = torch.vstack([X_dict[key], x_pad])

    return torch.stack(list(X_dict.values()))