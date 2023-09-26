import os
from mtz import MTZFile

TEST_DATA = os.path.join(os.path.dirname(__file__), "data")


def test_open():
    f = MTZFile(os.path.join(TEST_DATA, "scaled.mtz"))


def test_open_batches():
    f = MTZFile(os.path.join(TEST_DATA, "refl0.mtz"))
