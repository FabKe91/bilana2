import pytest

import numpy as np
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_equal,
    assert_array_equal,
)

import bilana2 as bl

from bilana2tests.datafiles import CHOL_INPUT

# =================================================================
#   This explicitly test attributes set for the first 5 frames of
#   dppc_chol20_300 simulation
# =================================================================

@pytest.fixture
def sysinfo():
    return bl.Systeminfo(CHOL_INPUT)

class TestInputAttributes:
    def test_system(self, sysinfo):
        assert sysinfo.system == "dppc_chol20", "System is set incorrectly"
    def test_temperature(self, sysinfo):
        assert sysinfo.temperature == 300, "Temperature is set incorrectly"
    def test_cutoff(self, sysinfo):
        assert sysinfo.cutoff == 1.0, "Cutoff is set incorrectly"
    def test_reference_atomselection(self, sysinfo):
        assert sysinfo.reference_atomselection ==\
            "(resname DPPC DUPC and name P) or (resname CHL1 and name O3) or (resname CHIM and name C20)",\
                "Reference string incorrectly read"
    def test_molecules(self, sysinfo):
        assert sysinfo.molecules == ["DPPC", "CHL1"], "Molecules are incorrect"
    #def test_mdfilespath(self, sysinfo):
    #    assert sysinfo.test_mdfilespath ==
    def test_forcefieldname(self, sysinfo):
        assert sysinfo.forcefieldname == "charmm", "Force field name set incorrectly"


class TestTimes:
    def test_t_start_corr(self, sysinfo):
        assert sysinfo.t_start == 0, "t start is set incorrectly"
    def test_t_end_corr(self, sysinfo):
        assert sysinfo.t_end == 5000, "t end is set incorrectly"
    def test_dt_corr(self, sysinfo):
        assert sysinfo.dt == 1000, "dt is set incorrectly"
    def test_t_start_inc(self, sysinfo):
        pass
    def test_t_end_inc(self, sysinfo):
        pass
    def test_dt_inc(self, sysinfo):
        pass
    def test_within_range(self, sysinfo):
        assert sysinfo.within_timerange(2000), "Is not found in range"
    def test_not_within_range_tmax(self, sysinfo):
        assert not sysinfo.within_timerange(1000000000000), "Is put in range though time>tmax"
    def test_not_within_range_dt(self, sysinfo):
        assert not sysinfo.within_timerange(500), "Is put in range though wrong dt"

class TestComposition:
    def test_all_lipids_found(self, sysinfo):
        assert sysinfo.number_of_lipids == 350, "Not all lipids found"
    def test_all_resids_set_correctly(self, sysinfo):
        assert np.all(np.array(sysinfo.lipid_resids) == np.arange(1, 351)), "Resids incorrect"
    def test_solvents_found(self, sysinfo):
        assert sysinfo.solvent_resnames == ["TIP3",], "not all solvent resnames found"
    def test_water_resname_correct(self, sysinfo):
        assert sysinfo.water_resname == "TIP3", "water set incorrectly"
    def test_pl_resnames_found(self, sysinfo):
        assert sysinfo.pl_resnames == ["DPPC",], "PL resname not found"
    def test_sterol_resnames_found(self, sysinfo):
        assert sysinfo.sterol_resnames == ["CHL1",], "PL resname not found"


class TestFilepaths:

    @pytest.fixture
    def path_class(self):
        pass

    def test_systempaths(self):
        pass
    def test_created_folders(self):
        pass


class TestConversions:

    @pytest.fixture
    def add_convert(self, sysinfo):
        return sysinfo.convert

    def test_resid_to_resname(self, add_convert):
        assert add_convert.resid_to_resname[1] == "DPPC", "DPPC not in resid_to_resname"
        assert add_convert.resid_to_resname[281] == "CHL1", "CHL1 not in resid_to_resname"

    def test_resid_to_leaflet(self, add_convert):
        assert add_convert.resid_to_leaflet[1] == 1, "resid 1 assigned incorrectly"
        assert add_convert.resid_to_leaflet[141] == 0, "leaflet 141 assigned incorrectly"
        assert add_convert.resid_to_leaflet[281] == 1, "leaflet 281 assigned incorrectly"
        assert add_convert.resid_to_leaflet[316] == 0 , "leaflet 316 assigned incorrectly"


# =================================================================
# =================================================================
# =================================================================