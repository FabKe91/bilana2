import pytest

import pandas as pd

import bilana2
from bilana2 import order

from ..datafiles import CHOL_INPUT, CHOL_ORDER, CHOL_NEIGHBOR, CHOL_TILT

CHOL_S_TEST         = "testfiles/scd_distribution_test.dat"
CHOL_SPROFILE_TEST  = "testfiles/scd_profile_test.dat"
CHOL_TILT_TEST      = "testfiles/tilt_test.csv"


@pytest.fixture
def get_sysinfo():
    return bilana2.Systeminfo(CHOL_INPUT)
@pytest.fixture
def get_neiblist():
    return bilana2.neighbor.get_neighbor_dict(neighborfilename=CHOL_NEIGHBOR)


class TestOrder:
    def test_order_parameter_is_correct(self, get_sysinfo, get_neiblist):
        order.create_cc_orderfiles(get_sysinfo, get_neiblist,
                outputfile_scd=CHOL_S_TEST,
                outputfile_s_profile=CHOL_SPROFILE_TEST,
                with_tilt_correction=CHOL_TILT,
                )
        dat_orig = pd.read_table(CHOL_ORDER, delim_whitespace=True)
        dat_test = pd.read_table(CHOL_S_TEST, delim_whitespace=True)
        dat_orig = dat_orig.drop(columns=["DPPC", "CHL1"])
        dat_test = dat_test.drop(columns=["DPPC", "CHL1"])
        assert dat_orig[dat_orig.Time <= 5000].equals(dat_test)

    def test_order_calculation(self):
        pass
    def test_order_profile(self):
        pass

def test_tilt_is_correct(get_sysinfo):
    order.calc_tilt(get_sysinfo, filename=CHOL_TILT_TEST)
    dat_orig = pd.read_csv(CHOL_TILT)
    dat_test = pd.read_csv(CHOL_TILT_TEST)
    assert dat_orig[dat_orig.time <= 5000].equals(dat_test)
