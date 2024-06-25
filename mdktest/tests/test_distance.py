import numpy as np
from numpy.testing import assert_allclose
import pytest

import MDAnalysis as mda
import MDAnalysisTests.datafiles as data

from mdktest import distance


# ADK FRET pairs
# see doi: 10.1016/j.jmb.2009.09.009.
fret_pairs = [
    ("resname ALA and resid 55", "resname VAL and resid 169"),
    ("resname ALA and resid 127", "resname ALA and resid 194"),
    ("resname ILE and resid 52", "resname LYS* and resid 145"),
]

@pytest.fixture
def u_nopbc():
    # simulation without periodic boundaries
    # (box == None)
    return mda.Universe(data.PSF, data.DCD)


@pytest.fixture
def u_pbc():
    # simulation with periodic boundaries and protein split
    # across unit cell
    return mda.Universe(data.TPR, data.XTC)


def _minimage_distance(u):
    # helper function
    a = sum(u.select_atoms(pair[0] + " and name CA") for pair in fret_pairs)
    b = sum(u.select_atoms(pair[1] + " and name CA") for pair in fret_pairs)

    r = distance.minimage_distance(a, b)

    return r

def test_minimage_distance_nopbc(u_nopbc):
    reference = np.array([12.91501288, 30.57278011, 28.92306815])
    r = _minimage_distance(u_nopbc)
    assert_allclose(r, reference)


def test_minimage_distance_pbc(u_pbc):
    reference = np.array([29.67992244, 38.36642607, 42.18877151])
    r = _minimage_distance(u_pbc)
    assert_allclose(r, reference)