
from .. import _utils as utils

import vinplots


class PlotUMAP(utils.ABCParse):
    def __init__(
        self,
        adata,
        title="10x PBMCs (~3k cells)",
        wspace=0.1,
        rm_ticks=True,
        spines_to_delete=["top", "bottom", "right", "left"],
        labels=["UMAP-1", "UMAP-2"],
        use_key="X_umap",
        cmap_key="cmap",
        alpha=0.65,
    ):
        self.__parse__(locals(), public=["adata"])

    def __config__(self):
        qiuck_plot_params = utils.FunctionInspector(func=vinplots.quick_plot)

        self.fig, self.axes = vinplots.quick_plot(
            nplots=1, ncols=1, **qiuck_plot_params(self._PARAMS)
        )

    @property
    def CMAP(self):
        return self.adata.uns[self._cmap_key]

    def _plot_grouped(self):
        for en, (group, group_df) in enumerate(self.adata.obs.groupby(self.groupby)):
            if group in self.plot_behind:
                z = 0
            else:
                z = en + 1
            xu = self.adata[group_df.index].obsm[self._use_key]
            self.axes[0].scatter(
                xu[:, 0],
                xu[:, 1],
                label=group,
                color=self.CMAP[group],
                alpha=self._alpha,
                zorder = z,
            )

    def _format(self):
        for ax in self.axes:
            ax.set_xlabel(self._labels[0])
            ax.set_ylabel(self._labels[1])
        self.axes[0].legend(loc=(1, 0.25), edgecolor="w")
        self.axes[0].set_title(self._title)

    def __call__(
        self,
        groupby=None,
        plot_behind = ["Undifferentiated"],
    ):
        self.__update__(locals())

        self.__config__()
        self._plot_grouped()
        self._format()