
__module_name__ = "__init__.py"
__doc__ = """Aux. module to organize AnnData/torch datasets into test/val/train subsets."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: ----------------------------------------------------------
import anndata
import inspect


# -- import local dependencies: ------------------------------------------------
from .._core._AnnDataset import AnnDataset
from ._split import split


# -- supporting functions: -----------------------------------------------------
def func_params(func):
    return list(inspect.signature(func).parameters.keys())


def extract_func_kwargs(func, kwargs):
    func_kwargs = {}
    params = func_params(func)
    for k, v in kwargs.items():
        if k in params:
            func_kwargs[k] = v
    return func_kwargs


# -- API-facing class: ---------------------------------------------------------
class AnnDatasetSplit:
    def __parse__(self, kwargs, ignore=["self"]):
        self.kwargs = {}
        for k, v in kwargs.items():
            if not k in ignore:
                setattr(self, k, v)
                self.kwargs[k] = v
                
    def __init__(
        self,
        adata: anndata.AnnData,
        use_key: str = "X_pca",
        groupby: str = None,
        obs_keys: str = None,
        percent_val: float = 0.2,
        train_key: str = "train",
        test_key: str = "test",
        **kwargs,
    ):
        self.__parse__(locals())
        self._obs_cols = adata.obs.columns.tolist()
        
    def _train_test_obs_keys(self):
        return all([key in self._obs_cols for key in [self.train_key, self.test_key]])

    def to_dataset(self, adata):
        
        AnnDataset_kwargs = extract_func_kwargs(AnnDataset, self.kwargs)
        
        return AnnDataset(
            adata,
            use_key=self.use_key,
            groupby=self.groupby,
            obs_keys=self.obs_keys,
            silent=True,
            **AnnDataset_kwargs
        )

    def on_test_train(self):
        """
        Split cells in adata on `train` / `test` columns in adata.obs.
        """
        if self._train_test_obs_keys():
            self.train_adata, self.test_adata = (
                self.adata[self.adata.obs[self.train_key]],
                self.adata[self.adata.obs[self.test_key]],
            )
            self.train_dataset = self.to_dataset(self.train_adata)
            self.test_dataset = self.to_dataset(self.test_adata)
        else:
            self.train_dataset, self.test_dataset = None, None

    def allocate_validation(self):
        """
        Split the previously allocated `train_dataset` into training and
        validation subsets.
        """
        if self.percent_val:
            self.train_dataset, self.val_dataset = split(
                self.train_dataset,
                percentages=[(1 - self.percent_val), self.percent_val],
            )
        else:
            self.val_dataset
