
__module_name__ = "_TimeResolvedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import local dependencies #
# ------------------------- #
from ._functions._use_X import _use_X
from ._AnnDataset import AnnDataset
from ._functions._pad_time_resolved_dataset import _pad_time_resolved_dataset


class TimeResolvedAnnDataset(AnnDataset):
    def __init__(self, adata, use_key="X_pca"):
        super(TimeResolvedAnnDataset, self).__init__(adata, use_key)

        self.X = _pad_time_resolved_dataset(self._adata)