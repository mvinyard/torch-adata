
__module_name__ = "_register_args.py"
__doc__ = """Argument registering module"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -----------------------------------------------------------------------
def register_attr_names(dataset, attr_names, obs_keys):

    """
    Notes:
    ------
    I also considered the option of using 'obs_{}'.format(N) as the
    key names. Instead, I decided to just replace and potential spaces
    in the obs_key passed with an underscore.
    """
    if not attr_names:
        attr_names = ["_".join(key.split(" ")) for key in obs_keys]
    if not isinstance(attr_names, list):
        attr_names = [attr_names]
    setattr(dataset, "_attr_names", ["X"] + attr_names)
    return attr_names


def register_one_hot_encode(
    one_hot: list([bool, ..., bool]),
    obs_keys: list([str, ..., str]),
)->list([bool, ..., bool]):
    """
    Parameters:
    -----------
    one_hot_encode
        a list, for each key in obs_keys, of boolean indicators that denote
        whether or not the corresponding column should be one-hot encoded.
        type: list

    obs_keys
        type: list([str, ..., str])

    Returns:
    --------
    one_hot_encode
        a list, for each key in obs_keys, of boolean indicators that denote
        whether or not the corresponding column should be one-hot encoded.
        type: list
    """
    if not one_hot:
        return [False for i in range(len(obs_keys))]
    if not isinstance(one_hot, list):
        return [one_hot]
    return one_hot


def register_args(dataset,
                  obs_keys: list([str, ..., str]),
                  attr_names: list([str, ..., str]),
                  one_hot: list([bool, ..., bool]),
                 ):
    """
    Parameters:
    -----------
    dataset

    obs_keys
        type: list([str, ..., str])

    attr_names
        type: list([str, ..., str])

    one_hot_encode
        a list, for each key in obs_keys, of boolean indicators that denote
        whether or not the corresponding column should be one-hot encoded.
        type: list([bool, ..., bool])

    Returns:
    --------
    attr_names
        type: list([str, ..., str])

    one_hot_encode
        a list, for each key in obs_keys, of boolean indicators that denote
        whether or not the corresponding column should be one-hot encoded.
        type: list([bool, ..., bool])
    """
    if obs_keys:
        attr_names = register_attr_names(dataset, attr_names, obs_keys)
        one_hot = register_one_hot_encode(one_hot, obs_keys)
        return attr_names, one_hot

    setattr(dataset, "_attr_names", ["X"])
    return None, None
