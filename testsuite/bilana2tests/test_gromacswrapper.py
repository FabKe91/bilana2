
import pytest

from bilana2 import _gmx

from bilana2tests.datafiles import CHOL_GRO, CHOL_TRJ

GMXCHECK_STDOUT = "Checking file /home/f_kell07/software/Bilana2/testsuite/bilana2tests/data/mdfiles_test/md_trj/dppc_chol20_300_whole.xtc\n"

def test_run_gmx():
    inp = ["-f", CHOL_TRJ]
    out, err = _gmx.exec_gromacs("gmx", "check", inp)
    assert out == GMXCHECK_STDOUT, "STDOUT is different"
def test_run_mdrun():
    pass
def test_err_code():
    pass


def test_log():
    pass
