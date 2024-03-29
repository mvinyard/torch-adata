
__module_name__ = "_group_and_pad_adata.py"
__doc__ = """Module for performing re-sampling upon grouping by a column in adata.obs"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import numpy as np
import anndata
import torch


NoneType = type(None)


def max_len(grouped_idx: dict) -> int:
    """
    Return the max len grouped index. Should match the dataset.__len__() attribute.
    Parameters:
    -----------
    grouped_idx
        type: dict
    Returns:
    --------
    max_len
        type: int
    """
    return max([v.shape[0] for v in grouped_idx.values()])


def _sample_indices_for_padding(
    grouped_idx: dict, sampling_replacement: bool = True, sampling_weights=None
) -> dict:

    """
    Parameters:
    -----------
    grouped_idx
        type: dict
    sampling_replacement
        type: bool
        default: True
    Returns:
    --------
    group_padding
        type: dict
    Notes:
    ------
    """

    group_padding = {}

    for group, group_idx in grouped_idx.items():
        n_pad = max_len(grouped_idx) - group_idx.shape[0]
        if n_pad:
            if not isinstance(sampling_weights, NoneType):
                p_ = sampling_weights[group_idx]
                p = p_ / p_.sum()
            else:
                p = None
            group_padding[group] = np.random.choice(
                group_idx, size=n_pad, replace=sampling_replacement, p=p
            )

    return group_padding


def _pad_indices(grouped_idx: dict, sampling_replacement=True, sampling_weights=None) -> dict:

    """
    Parameters:
    -----------
    grouped_idx
        type: dict
    Returns:
    --------
    grouped_idx
        type: dict
    Notes:
    ------
    """

    padding_idx = _sample_indices_for_padding(
        grouped_idx,
        sampling_replacement,
        sampling_weights,
    )

    for group, group_pad_idx in padding_idx.items():
        grouped_idx[group] = np.concatenate([grouped_idx[group], group_pad_idx])

    return grouped_idx


def group_and_pad_adata(adata: anndata.AnnData, groupby: str, sampling_replacement=True, sampling_weights=None) -> dict:

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData
    groupby
        type: str
    Returns:
    --------
    grouped_idx
    type: dict
    Notes:
    ------
    """

    grouped = adata.obs.groupby(groupby)
    idx = _pad_indices(grouped.indices, sampling_replacement=True, sampling_weights=sampling_weights)
    groups = torch.Tensor(list(grouped.groups.keys()))

    return idx, groups
