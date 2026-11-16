"""Shared observation vocabularies, established before slicing or sampling."""

import pandas as pd
from sklearn.preprocessing import OneHotEncoder
import torch

from ._data_typing import as_list


def categorical_observation(values, one_hot=False):
    # Numeric pandas categoricals historically travel as numeric observations.
    dtype = values.to_numpy().dtype if isinstance(values.dtype, pd.CategoricalDtype) else values.dtype
    return (one_hot or pd.api.types.is_object_dtype(dtype)
            or pd.api.types.is_string_dtype(dtype))


def category_codes(values, categories, key):
    if values.isna().any():
        raise ValueError(f"Observation {key!r} contains missing categorical labels")
    codes = pd.Index(categories).get_indexer(values.to_numpy())
    if (codes < 0).any():
        unknown = values.iloc[codes < 0].unique().tolist()
        raise ValueError(f"Observation {key!r} contains unknown categorical labels: {unknown!r}")
    return codes


def build_category_vocabularies(adata, obs_keys, one_hot=False, category_vocabularies=None):
    """Copy supplied column orders and infer missing vocabularies from all rows.

    Inferred categories are sorted, observed values only. Missing categorical
    labels and values outside a supplied vocabulary are rejected.
    """
    keys = as_list(obs_keys)
    flags = one_hot if isinstance(one_hot, list) else [one_hot] * len(keys)
    if len(flags) != len(keys):
        raise ValueError("one_hot must contain one flag per observation key")
    vocabularies = {}
    for key, categories in (category_vocabularies or {}).items():
        categories = tuple(categories)
        index = pd.Index(categories)
        if not len(index) or index.hasnans or not index.is_unique:
            raise ValueError(f"Vocabulary for {key!r} must contain unique, non-missing categories")
        vocabularies[key] = categories
    for key, flag in zip(keys, flags):
        values = adata.obs[key]
        if key not in vocabularies and categorical_observation(values, flag):
            if values.isna().any():
                raise ValueError(f"Observation {key!r} contains missing categorical labels")
            encoder = OneHotEncoder().fit(values.to_numpy().reshape(-1, 1))
            vocabularies[key] = tuple(encoder.categories_[0].tolist())
        if key in vocabularies:
            category_codes(values, vocabularies[key], key)
    return vocabularies


def encode_categories(values, categories, key):
    codes = category_codes(values, categories, key)
    return torch.nn.functional.one_hot(
        torch.as_tensor(codes, dtype=torch.long), num_classes=len(categories)
    ).to(torch.float32)
