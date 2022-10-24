
__module_name__ = "_return_data_on_axis.py"
__doc__ = """Key method for returning stacked data given 'groupby' (or not) arguments."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify package version: ------------------------------------------------------------
__version__ = "0.0.17"


# ----------------------------------------------------------------------------------------
def return_data_on_axis(dataset, idx):
    """
    Automatically return fetched / sampled batches of data along the right data axes.
    """
    if dataset._data_axis:
        return [dataset.groups.unique()] + [
            getattr(dataset, key)[:, idx] for key in dataset._attr_names
        ]
    return [getattr(dataset, key)[idx] for key in dataset._attr_names]