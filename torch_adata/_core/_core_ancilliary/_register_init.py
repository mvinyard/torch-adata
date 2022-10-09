
__module_name__ = "_register_init.py"
__doc__ = """
          Module for registering all arguments as well as parsing data during the init setup.
          """
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
from ._fetch_data import Fetch
from ._identity_msg import identity_msg


# -- supporting functions: ---------------------------------------------------------------
def rm_space(key, filler="_"):
    return filler.join(key.split(" "))

def as_list(item):
    if item:
        if not isinstance(item, list):
            return [item]
        return item
    return []

def compile_attributes(obs_keys, aux_keys):
    obs_keys = as_list(obs_keys)
    aux_keys = as_list(aux_keys)
    attr_keys = as_list(obs_keys) + as_list(aux_keys)
    return attr_keys, obs_keys, aux_keys

def count_attrs(dataset):
    setattr(dataset, "_n_attrs", len(dataset._attr_names))
    
def register_attr_names(dataset, attr_names, given_keys):

    """
    Notes:
    ------
    I also considered the option of using 'obs_{}'.format(N) as the
    key names. Instead, I decided to just replace and potential spaces
    in the obs_key passed with an underscore.
    """

    if not attr_names:
        attr_names = [rm_space(key) for key in given_keys]
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
                  aux_keys: list([str, ..., str]),
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
        
    Notes:
    ------
    (1) obs and aux attribute keys are combined into a single list. obs_keys
        listed before aux_keys
    """

    
    attr_keys, obs_keys, aux_keys = compile_attributes(obs_keys, aux_keys)

    if attr_keys:
        attr_names = register_attr_names(dataset, attr_names, attr_keys)
        one_hot = register_one_hot_encode(one_hot, obs_keys)
        print(" - ONE HOT: ", one_hot)
        return attr_keys, obs_keys, aux_keys, attr_names, one_hot

    setattr(dataset, "_attr_names", ["X"])
    return [None] * 3


# -- Main module function: ---------------------------------------------------------------
def register_init(
    dataset,
    adata,
    groupby,
    use_key,
    obs_keys,
    aux_keys,
    attr_names,
    one_hot,
    silent,
):
    """
    Register passed inputs to AnnDataset class.

    Notes:
    ------
    As a rule, we'll always keep `X` as a required attribute. Optional attributes then
    come in two additional forms: (1) `X_like` and (2) obs with the following shapes:
    (1) [n_groups, n_samples, n_dim] (same as `X`)
    (2) [n_groups, n_samples, 1]
    """

    dataset._silent = silent

    attr_keys, obs_keys, aux_keys, attr_names, one_hot = register_args(
        dataset, obs_keys, aux_keys, attr_names, one_hot
    )

    fetch = Fetch(adata)

    if groupby:
        dataset._data_axis = 1
        dataset._grouped_by = groupby
        dataset.X, obs_data, aux_data = fetch.grouped_adata(
            groupby=groupby,
            use_key=use_key,
            obs_keys=obs_keys,
            aux_keys=aux_keys,
            attr_names=attr_names,
            one_hot=one_hot,
        )
        fetch.update_attrs(dataset, obs_data, aux_data)

    else:
        dataset._data_axis = 0
        dataset._grouped_by = None
        dataset.X = fetch.X(use_key)
        if obs_keys:
            obs_data = fetch.multi_obs(obs_keys, attr_names, one_hot)
            fetch.update_obs_attrs(dataset, obs_data)
        if len(attr_keys) > len(obs_keys): # if True, has aux attrs
            aux_data = {}
            for k in range(len(obs_keys), len(attr_keys)):
                aux_attr = attr_names[k]
                print(" - AUX: {}".format(aux_attr))
                setattr(dataset, aux_attr, fetch.X(aux_attr))
#                 aux_data[aux_attr] = fetch.X(aux_attr)
        
    count_attrs(dataset)
    if not dataset._silent:
        print(identity_msg(dataset))