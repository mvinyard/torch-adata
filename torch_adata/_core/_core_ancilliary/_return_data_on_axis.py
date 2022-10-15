
def return_data_on_axis(dataset, idx):
    """
    Automatically return fetched / sampled batches of data along the right data axes.
    """
    if dataset._data_axis:
        return [dataset.groups.unique()] + [
            getattr(dataset, key)[:, idx] for key in dataset._attr_names
        ]
    return [getattr(dataset, key)[idx] for key in dataset._attr_names]