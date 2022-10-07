
__module_name__ = "_idx.py"
__doc__ = """Utility module to generate an index within the bounds of len(dataset)"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import numpy as np


# -- sole function: ----------------------------------------------------------------------
def idx(
    dataset,
    size: int = None,
    replace: bool = False,
) -> np.ndarray:

    """
    Create a an index, sampling from range(len(dataset)).

    Parameters:
    -----------
    dataset

    size
        type: int

    replace
        type: bool

    Returns:
    --------
    idx
        type: numpy.ndarray
    """

    if not size:
        size = int(len(dataset) / 100)
    return np.random.choice(range(len(dataset)), size=size, replace=replace)
