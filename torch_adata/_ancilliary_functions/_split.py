
__module_name__ = "_split.py"
__doc__ = """Helper function(s) for splitting the dataset."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import torch


# -- supporting functions: ---------------------------------------------------------------
def uniform_group_sizes(
    dataset: torch.utils.data.Dataset, n_groups: int = 2
) -> list([int, ..., int]):
    """Split dataset length with even proportions according to number of groups"""

    div, mod = divmod(len(dataset), n_groups)
    return [div] * (n_groups - 1) + [div + mod]


def proportioned_group_sizes(
    dataset: torch.utils.data.Dataset,
    percentages: list([float, ..., float]),
    remainder_idx: int = -1,
) -> list([int, ..., int]):
    """Split dataset length with specified proportions/number of groups"""
    n_groups = len(percentages)
    if sum(percentages) < 1:
        percentages.append(round(1 - sum(percentages), 5))

    split_lengths = [int(len(dataset) * pct_i) for pct_i in percentages]
    remainder = len(dataset) - sum(split_lengths)
    split_lengths[remainder_idx] += remainder

    return split_lengths


def calculate_split_lengths(
    dataset: torch.utils.data.Dataset,
    n_groups: int = 2,
    percentages: list([float, ..., float]) = None,
) -> list([int, ..., int]):
    """Split the length of the dataset into proportioned subsets"""
    if percentages:
        return proportioned_group_sizes(dataset, percentages, remainder_idx=-1)
    return uniform_group_sizes(dataset, n_groups)


# -- main module function: ---------------------------------------------------------------
def split(
    dataset: torch.utils.data.Dataset,
    n_groups: int = 2,
    percentages: list([float, ..., float]) = None,
) -> list([torch.utils.data.Dataset, ..., torch.utils.data.Dataset]):
    """
    Split dataset using torch.utils.data.random_split.

    Parameters:
    -----------
    dataset
        type: torch.utils.data.Dataset

    n_groups
        type: int
        default: 2

    percentages
        type: list([float, ..., float])
        default: None

    Returns:
    --------
    list([torch.utils.data.Dataset, ..., torch.utils.data.Dataset])

    Notes:
    ------
    (1) Uses the torch.utils.data.random_split function to actually do the split.
    """
    split_lengths = calculate_split_lengths(dataset, n_groups, percentages)
    return torch.utils.data.random_split(dataset, lengths=split_lengths)