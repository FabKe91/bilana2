import pytest

import bilana2

import pandas as pd

from bilana2tests.datafiles import CHOL_INPUT, CHOL_RESINDEX, CHOL_NEIGHBOR

CHOL_RESINDEX_TEST  = "testfiles/resindex_all_test.ndx"
CHOL_ALLENERGY_TEST = "testfiles/all_energies_test.ndx"

@pytest.fixture
def get_sysinfo():
    return bilana2.Systeminfo(CHOL_INPUT)
@pytest.fixture
def get_neiblist():
    return bilana2.neighbor.get_neighbor_dict(neighborfilename=CHOL_NEIGHBOR)
@pytest.fixture
def get_energy(get_neiblist):
    energy = bilana2.Energy("complete", get_neiblist,
            inputfilepath=CHOL_INPUT,
            resindex_all=CHOL_RESINDEX)
    return energy

def test_loadenergy(get_neiblist):
    energy = bilana2.Energy("complete", get_neiblist,
            inputfilepath=CHOL_INPUT,
            resindex_all=CHOL_RESINDEX)

class TestEnergyCalculation:

    @pytest.mark.skip(reason="lazy")
    def test_energy_res1(self, get_energy):
        get_energy.run_calculation(resids=[1])
        bilana2.energy.write_energyfile(get_energy)
        dat_orig = pd.read_csv(get_energy.all_energies)
        dat_test = pd.read_csv(CHOL_ALLENERGY_TEST)
        assert dat_orig.equals(dat_test)


    def test_mdp_res1(self, get_energy):
        pass
    def test_edr_res1(self, get_energy):
        pass
    def test_xvg_res1(self, get_energy):
        pass
    def test_all_energies(self, get_energy):
        pass

class TestIndexCreation:
    @pytest.mark.skip(reason="found no easy way of testing it")
    def test_resindex_all(self, get_sysinfo):
        bilana2.energy.create_indexfile(get_sysinfo,
            resindex_filename=CHOL_RESINDEX_TEST)
        dat_orig = pd.read_csv(CHOL_RESINDEX)
        dat_test = pd.read_csv(CHOL_RESINDEX_TEST)
        assert dat_orig.equals(dat_test)
