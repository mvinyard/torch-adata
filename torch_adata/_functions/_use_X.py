
__module_name__ = "_use_X.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import torch


def _use_X(adata, use_key="X_pca"):

    if use_key == "X":
        return torch.Tensor(adata.X.toarray())

    elif use_key in adata.obsm_keys():
        return torch.Tensor(adata.obsm[use_key])

    else:
        print("No suitable array found for dataset...")