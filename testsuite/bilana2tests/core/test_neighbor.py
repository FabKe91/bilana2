import pytest

import pandas as pd

import bilana2
from bilana2 import neighbor

from ..datafiles import CHOL_INPUT, CHOL_NEIGHBOR

NEIGHBOR_INFO_TEST = "testfiles/neighbor_info_test"

@pytest.fixture
def get_sysinfo():
    return bilana2.Systeminfo(CHOL_INPUT)

def test_get_neighbor_dict():
    neib_dict = neighbor.get_neighbor_dict(neighborfilename=CHOL_NEIGHBOR)
    assert neib_dict[0][1] == [2, 3, 5, 7], "Loading the dict failed"

def test_get_ref_positions():
    pass

@pytest.mark.skip(reason="new neighbor_info routine will bring different results")
def test_neighbor_info_correct(get_sysinfo):
    neighbor.write_neighbor_info(get_sysinfo, outputfilename=NEIGHBOR_INFO_TEST)
    dat_orig = pd.read_table(CHOL_NEIGHBOR, delim_whitespace=True)
    dat_test = pd.read_table(NEIGHBOR_INFO_TEST, delim_whitespace=True)
    assert dat_orig[dat_orig.Time < 6000].equals(dat_test.sort_values(by=["Time", "Resid"])),\
            "new neighbor file is different"
