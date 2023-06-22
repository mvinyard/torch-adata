
import vinplots
import anndata
import os
import web_files

from .. import _utils as utils
from ._plot_umap import PlotUMAP

class LARRY(utils.ABCParse):
    _HTTP_ADDRESS = "https://figshare.com/ndownloader/files/38171943"

    def __init__(self, filename="adata.LARRY.h5ad", data_dir="data", silent=False):
        self._config(locals())

    def _config(self, kwargs):
        """"""

        self.__parse__(kwargs, public=[None])

        self._INFO = utils.InfoMessage()

        if not os.path.exists(self._data_dir):
            os.mkdir(self._data_dir)

    def _download(self):
        msg = "Downloading {:.<10}...".format(self._filename)
        self._INFO(msg, end=" ")
        f = web_files.WebFile(
            http_address=self._HTTP_ADDRESS,
            local_path=self.local_path,
        )
        f.download()
        print("Done.")

    @property
    def local_path(self):
        return os.path.join(self._data_dir, self._filename)

    @property
    def adata(self):
        if not hasattr(self, "_adata"):
            if not os.path.exists(self.local_path):
                self._download()
            self._adata = anndata.read_h5ad(self.local_path)
            self._adata.uns["cmap"] = vinplots.colors.LARRY_in_vitro
            if not self._silent:
                print(self._adata)
                self.plot()
        return self._adata
    
    def plot(self):
        umap_plot = PlotUMAP(
            self._adata, title="LARRY hematopoiesis"
        )
        umap_plot(groupby="Cell type annotation")
