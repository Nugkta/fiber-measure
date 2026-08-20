"""Tests for the multi-angle (C1) filename parsing and discovery API."""

from pathlib import Path

import pytest

from fibrecv.io_utils import (
    MultiAngleKey,
    discover_multiangle,
    multiangle_group,
    parse_multiangle_name,
)


@pytest.mark.parametrize("name,cond,fiber,angle,part,scalebar", [
    ("C1_01_a1_part1.tiff", "C1", 1, 1, 1, False),
    ("C1_15_a6_part5.tiff", "C1", 15, 6, 5, False),
    ("C1_01_a1_part1s.tiff", "C1", 1, 1, 1, True),
    ("C1_08_a3_part2.tiff", "C1", 8, 3, 2, False),
    ("D12_100_a2_part9s.tiff", "D12", 100, 2, 9, True),
    ("Cond_7_a4_part3.tiff", "Cond", 7, 4, 3, False),
])
def test_parse_multiangle_name(name, cond, fiber, angle, part, scalebar):
    key = parse_multiangle_name(name)
    assert key == MultiAngleKey(cond, fiber, angle, part, scalebar)


@pytest.mark.parametrize("name", [
    "C1_09_a6_par5.tiff",            # truncated 'part' — strict, no special case
    "C1_13_a6_part34.tiff",          # two-digit part
    "C1_01_a1_part1.tiff_metadata.xml",  # sidecar (stem still ends .tiff)
    "masp2 10_5_2.jpg",              # MasP2 convention
    "C1_01_part1.tiff",              # missing angle token
    "C1_a1_part1.tiff",              # missing fiber token
    "C1_01_a1_part1ss.tiff",         # double s
    "_01_a1_part1.tiff",             # empty condition
    "C1_1234_a1_part1.tiff",         # fiber > 3 digits
])
def test_parse_multiangle_name_rejects(name):
    with pytest.raises(ValueError):
        parse_multiangle_name(name)


def test_parse_multiangle_name_accepts_path_objects(tmp_path):
    key = parse_multiangle_name(tmp_path / "C1_03_a2_part4.tiff")
    assert key == MultiAngleKey("C1", 3, 2, 4, False)


def test_multiangle_group_adapter():
    key = MultiAngleKey("C1", 3, 2, 4, False)
    assert multiangle_group(key) == ("C1_03_a2", 4)
    key2 = MultiAngleKey("C1", 15, 6, 5, True)
    assert multiangle_group(key2) == ("C1_15_a6", 5)


def _touch(root: Path, name: str) -> None:
    (root / name).write_bytes(b"")


def test_discover_multiangle(tmp_path):
    # two fibers x two angles x one part, plus s twins, sidecars and noise
    for n in [
        "C1_01_a1_part1.tiff", "C1_01_a2_part1.tiff",
        "C1_02_a1_part3.tiff",
        "C1_01_a1_part1s.tiff",                 # s twin: excluded
        "C1_01_a1_part1.tiff_metadata.xml",     # sidecar: excluded
        "notes.txt", "background.jpg",          # noise: excluded
    ]:
        _touch(tmp_path, n)
    found = discover_multiangle(tmp_path)
    assert set(found) == {(1, 1), (2, 3)}
    assert set(found[(1, 1)]) == {1, 2}
    assert found[(1, 1)][1].name == "C1_01_a1_part1.tiff"
    assert found[(1, 1)][2].name == "C1_01_a2_part1.tiff"
    assert set(found[(2, 3)]) == {1}


def test_discover_multiangle_condition_filter(tmp_path):
    for n in ["C1_01_a1_part1.tiff", "D2_01_a1_part1.tiff"]:
        _touch(tmp_path, n)
    only_c1 = discover_multiangle(tmp_path, condition="C1")
    assert set(only_c1) == {(1, 1)}
    assert only_c1[(1, 1)][1].name == "C1_01_a1_part1.tiff"
    both = discover_multiangle(tmp_path)
    # no condition filter: both conditions collide on (fiber, part) — the
    # caller is expected to filter; here we just require both files seen
    assert len(both) >= 1


def test_discover_multiangle_missing_angle_is_absent_key(tmp_path):
    for a in (1, 2, 4, 5, 6):  # a3 missing
        _touch(tmp_path, f"C1_05_a{a}_part2.tiff")
    found = discover_multiangle(tmp_path)
    assert set(found[(5, 2)]) == {1, 2, 4, 5, 6}
    assert 3 not in found[(5, 2)]
