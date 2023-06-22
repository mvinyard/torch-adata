
# -- import packages: --------------------------------------------------------
import scanpy as sc
import pandas as pd
import vinplots


from .. import _utils as utils
from ._plot_umap import PlotUMAP

class PBMC3k(utils.ABCParse):
    def __init__(self, silent=False):
        self.__parse__(locals(), public=[None])

    def _configure(self):
        _adata = sc.datasets.pbmc3k_processed()
        cell_type_df = pd.DataFrame(_adata.obs["louvain"].unique()).reset_index()
        cell_type_df.columns = ["cell_type_idx", "louvain"]
        df_obs = _adata.obs.merge(cell_type_df, on="louvain", how="left")
        df_obs.index = df_obs.index.astype(str)
        _adata.obs = df_obs
        _adata.uns["cmap"] = vinplots.colors.pbmc3k
        return _adata

    @property
    def adata(self):
        if not hasattr(self, "_adata"):
            self._adata = self._configure()
            if not self._silent:
                print(self._adata)
                self.plot()
        return self._adata

    def plot(self):
        umap_plot = PlotUMAP(
            self._adata, title="10x PBMCs (~3k cells)"
        )
        umap_plot(groupby="louvain")
        
