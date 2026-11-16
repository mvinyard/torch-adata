import anndata
import numpy as np
import pandas as pd
import pytest
import torch

from torch_adata import AnnDataset, LightningAnnDataModule


def make_adata(labels, dtype="object", times=None):
    n = len(labels)
    obs = pd.DataFrame({
        "t": times if times is not None else [0] * (n // 2) + [1] * (n - n // 2),
        "label": pd.Series(labels, dtype=dtype),
        "other": ["x", "y"] * (n // 2),
        "weight": np.arange(n, dtype=float) + 0.5,
    })
    obs.index = obs.index.astype(str)
    return anndata.AnnData(np.arange(n, dtype=np.float32).reshape(-1, 1), obs=obs)


@pytest.mark.parametrize("dtype", ["object", "category"])
@pytest.mark.parametrize("labels", [
    ["A", "B", "A", "B", "A", "B", "C", "A"],
    ["A", "B", "A", "B", "B", "C", "B", "C"],
])
def test_group_columns_have_global_meaning(dtype, labels):
    adata = make_adata(labels, dtype)
    ds = AnnDataset(adata, "X", groupby="t", obs_keys=["label", "other", "weight"], silent=True)
    expected = {"A": [1, 0, 0], "B": [0, 1, 0], "C": [0, 0, 1]}
    assert ds.label.shape == (2, 4, 3)
    assert ds.category_vocabularies["label"] == ("A", "B", "C")
    for group in range(2):
        for row in range(4):
            cell = int(ds.X[group, row, 0])
            assert ds.label[group, row].tolist() == expected[labels[cell]]
            assert ds.other[group, row].tolist() == ([1, 0] if cell % 2 == 0 else [0, 1])
            assert ds.weight[group, row].tolist() == [cell + 0.5]


@pytest.mark.parametrize("dtype", ["object", "category"])
def test_lightning_shares_schema_without_changing_split_membership(dtype):
    adata = make_adata(["A", "B", "B", "C", "C", "C"], dtype, [0, 1, 0, 1, 0, 1])
    for name, ids in {"train": [0, 1], "val": [2, 3], "test": [4, 5]}.items():
        adata.obs[name] = np.isin(np.arange(6), ids)
    original = adata.obs.copy(deep=True)
    dm = LightningAnnDataModule(adata, use_key="X", groupby="t", obs_keys=["label"], num_workers=0)
    expected = {"train": [[1, 0, 0], [0, 1, 0]], "val": [[0, 1, 0], [0, 0, 1]], "test": [[0, 0, 1], [0, 0, 1]]}
    for name, ids in {"train": [0, 1], "val": [2, 3], "test": [4, 5]}.items():
        ds = getattr(dm, name + "_dataset")
        assert ds.X.flatten().tolist() == ids
        assert ds.label[:, 0].tolist() == expected[name]
    pd.testing.assert_frame_equal(adata.obs, original)


def test_explicit_schema_reused_for_ungrouped_subset():
    adata = make_adata(["A", "B", "B", "C"])
    full = AnnDataset(adata, "X", obs_keys=["label"], silent=True)
    subset = AnnDataset(adata[[2, 3]], "X", obs_keys=["label"], silent=True,
                        category_vocabularies=full.category_vocabularies)
    assert subset.label.tolist() == [[0, 1, 0], [0, 0, 1]]


def test_vocabulary_exists_before_padding():
    adata = make_adata(["A", "B", "C", "C", "C", "C"], times=[0, 0, 1, 1, 1, 1])
    # Zero probability prevents B from appearing in the added rows.
    adata.obs["sampling"] = [1., 0., 1., 1., 1., 1.]
    ds = AnnDataset(adata, "X", groupby="t", obs_keys=["label"],
                    sampling_weight_key="sampling", silent=True)
    assert ds.label.shape == (2, 4, 3)
    assert ds.label[0].tolist() == [[1, 0, 0], [0, 1, 0], [1, 0, 0], [1, 0, 0]]
    assert ds.label[1].tolist() == [[0, 0, 1]] * 4


@pytest.mark.parametrize("missing", [None, np.nan, pd.NA])
def test_missing_categorical_labels_raise(missing):
    adata = make_adata(["A", missing, "A", "B"])
    with pytest.raises(ValueError, match="label.*missing"):
        AnnDataset(adata, "X", obs_keys=["label"], silent=True)


def test_unknown_label_raises_instead_of_remapping():
    adata = make_adata(["A", "B", "C", "C"])
    with pytest.raises(ValueError, match="label.*unknown"):
        AnnDataset(adata, "X", obs_keys=["label"], silent=True,
                   category_vocabularies={"label": ("A", "B")})


def test_explicit_numeric_one_hot_uses_global_categories():
    adata = make_adata(["A", "B", "C", "C"])
    adata.obs["code"] = [1, 2, 2, 3]
    ds = AnnDataset(adata, "X", groupby="t", obs_keys=["code", "weight"],
                    one_hot=[True, False], silent=True)
    assert ds.code.tolist() == [[[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [0, 0, 1]]]
    torch.testing.assert_close(ds.weight, torch.tensor([[[0.5], [1.5]], [[2.5], [3.5]]]))


@pytest.mark.parametrize("categories", [[], ["A", "A"], ["A", None]])
def test_invalid_vocabularies_raise(categories):
    with pytest.raises(ValueError, match="Vocabulary.*label"):
        AnnDataset(make_adata(["A", "B", "A", "B"]), "X", obs_keys=["label"],
                   category_vocabularies={"label": categories}, silent=True)


def test_supplied_order_is_preserved_and_copied():
    schema = {"label": ["C", "B", "A"]}
    ds = AnnDataset(make_adata(["A", "B", "A", "B"]), "X", obs_keys=["label"],
                    category_vocabularies=schema, silent=True)
    schema["label"].reverse()
    assert ds.category_vocabularies["label"] == ("C", "B", "A")
    assert ds.label.tolist() == [[0, 0, 1], [0, 1, 0], [0, 0, 1], [0, 1, 0]]


def test_lightning_explicit_schema_and_repeated_construction():
    adata = make_adata(["A", "B", "C", "C"])
    adata.obs["train"] = [True, True, False, False]
    dm = LightningAnnDataModule(adata, use_key="X", obs_keys=["label"], num_workers=0,
                                category_vocabularies={"label": ["C", "B", "A"]})
    for _ in range(2):
        assert dm.train_dataset.label.tolist() == [[0, 0, 1], [0, 1, 0]]
        assert dm.predict_dataset.label.tolist() == [[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 0, 0]]


def test_declared_unused_categories_do_not_change_inferred_columns():
    adata = make_adata(["A", "B", "A", "B"], "category")
    adata.obs["label"] = adata.obs["label"].cat.set_categories(["C", "B", "A"], ordered=True)
    ds = AnnDataset(adata, "X", obs_keys=["label"], silent=True)
    assert ds.label.tolist() == [[1, 0], [0, 1], [1, 0], [0, 1]]


def test_boolean_one_hot_applies_to_all_requested_columns():
    adata = make_adata(["A", "B", "A", "B"])
    adata.obs["code"] = [1, 2, 1, 2]
    ds = AnnDataset(adata, "X", obs_keys=["code", "label"], one_hot=True, silent=True)
    assert ds.code.tolist() == ds.label.tolist() == [[1, 0], [0, 1], [1, 0], [0, 1]]


def test_empty_subset_retains_explicit_width():
    adata = make_adata(["A", "B", "A", "B"])
    ds = AnnDataset(adata[:0], "X", obs_keys=["label"],
                    category_vocabularies={"label": ["A", "B"]}, silent=True)
    assert ds.label.shape == (0, 2)


def test_numeric_pandas_categories_retain_numeric_behavior():
    adata = make_adata(["A", "B", "A", "B"])
    adata.obs["code"] = pd.Categorical([1, 2, 2, 3])
    ds = AnnDataset(adata, "X", groupby="t", obs_keys=["code"], silent=True)
    assert ds.code.tolist() == [[[1], [2]], [[2], [3]]]
