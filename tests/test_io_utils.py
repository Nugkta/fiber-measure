"""Tests for filename parsing and image discovery (io_utils)."""

import pytest

from fibrecv.io_utils import discover_images, natural_key, parse_name


@pytest.mark.parametrize("name,group,rep", [
    ("masp2 10_5_2.jpg", "10_5", 2),           # current convention, unchanged
    ("masp2 3_1_10.jpg", "3_1", 10),           # multi-digit replicate
    ("MASP2 3_1_2.jpg", "3_1", 2),             # prefix case is irrelevant
    ("3-1-2.png", "3_1", 2),                   # dash separators, no prefix
    ("sampleA 10_5_2.tif", "10_5", 2),         # arbitrary text prefix
    ("fiber_3_1 (2).jpg", "3_1", 2),           # parenthesised replicate
    ("3_1_2.jpg", "3_1", 2),                   # bare numbers
    ("IMG_0123.jpg", "IMG", 123),              # single trailing number
    ("scan 7.jpeg", "scan", 7),                # single number, space separator
])
def test_parse_name(name, group, rep):
    assert parse_name(name) == (group, rep)


@pytest.mark.parametrize("name", [
    "background.jpg",        # no digits at all
    "masp2 3_1_2 copy.jpg",  # text after the trailing numbers
    "3.jpg",                 # single number with no prefix -> no group
])
def test_parse_name_rejects(name):
    with pytest.raises(ValueError):
        parse_name(name)


def test_natural_key_orders_numeric_groups():
    groups = ["10_5", "3_3", "3_1", "IMG"]
    assert sorted(groups, key=natural_key) == ["3_1", "3_3", "10_5", "IMG"]


def test_discover_images_filters_extensions_and_sidecars(tmp_path):
    names = [
        "masp2 3_1_2.jpg", "3-1-2.PNG", "scan 7.tif", "IMG_0123.jpeg",
        "notes.txt", "masp2 3_1_2.jpg_metadata.xml",
    ]
    for n in names:
        (tmp_path / n).write_bytes(b"")
    found = {p.name for p in discover_images(tmp_path)}
    assert found == {"masp2 3_1_2.jpg", "3-1-2.PNG", "scan 7.tif",
                     "IMG_0123.jpeg"}


def test_discover_images_natural_order(tmp_path):
    for n in ["masp2 10_1_1.jpg", "masp2 3_1_10.jpg", "masp2 3_1_2.jpg"]:
        (tmp_path / n).write_bytes(b"")
    assert [p.name for p in discover_images(tmp_path)] == [
        "masp2 3_1_2.jpg", "masp2 3_1_10.jpg", "masp2 10_1_1.jpg"]


def test_select_images_groups_filter(tmp_path):
    import argparse

    from fibrecv.run_measure import select_images

    for n in ["masp2 3_1_1.jpg", "masp2 3_1_2.jpg", "sampleA 10_5_1.png",
              "background.jpg"]:
        (tmp_path / n).write_bytes(b"")
    args = argparse.Namespace(root=str(tmp_path), glob=None, groups=["10_5"])
    assert [p.name for p in select_images(args)] == ["sampleA 10_5_1.png"]
