import MDAnalysis as mda
from mdktest import info
from MDAnalysisTests import datafiles as data

import pytest

@pytest.fixture
def u_charmm():
    udata = {
        "n_atoms": 3341,
        "Lx": 0,
        "Ly": 0,
        "Lz": 0,
        "alpha": 0,
        "beta": 0,
        "gamma": 0,
        "n_frames": 98,
        "totaltime": 96.9999914562418,
        "dt": 0.9999999119200186,
    }
    return mda.Universe(data.PSF, data.DCD), udata


@pytest.fixture
def u_pdb():
    udata = {
        "n_atoms": 47681,
        "Lx": 80.017,
        "Ly": 80.017,
        "Lz": 80.017,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 90.0,
        "n_frames": 1,
        "totaltime": 0.0,
        "dt": 1.0,
    }
    return mda.Universe(data.PDB), udata

@pytest.fixture(scope='module')
def u_gmx():
    reference = {
        "n_atoms": 47681,
        "Lx": 80.017006,
        "Ly": 80.017006,
        "Lz": 80.017006,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 90.0,
        "n_frames": 10,
        "totaltime": 900.0000686645508,
        "dt": 100.00000762939453,
    }
    universe = mda.Universe(data.TPR, data.XTC)
    return universe, reference

@pytest.mark.parametrize(
    "topology,trajectory,reference",
    [
        (
            data.TPR,
            data.XTC,
            {
                "n_atoms": 47681,
                "Lx": 80.017006,
                "Ly": 80.017006,
                "Lz": 80.017006,
                "alpha": 60.0,
                "beta": 60.0,
                "gamma": 90.0,
                "n_frames": 10,
                "totaltime": 900.0000686645508,
                "dt": 100.00000762939453,
            },
        ),
        (
            data.PSF,
            data.DCD,
            {
                "n_atoms": 3341,
                "Lx": 0,
                "Ly": 0,
                "Lz": 0,
                "alpha": 0,
                "beta": 0,
                "gamma": 0,
                "n_frames": 98,
                "totaltime": 96.9999914562418,
                "dt": 0.9999999119200186,
            },
        ),
        (
            data.PDB,
            None,
            {
                "n_atoms": 47681,
                "Lx": 80.017,
                "Ly": 80.017,
                "Lz": 80.017,
                "alpha": 60.0,
                "beta": 60.0,
                "gamma": 90.0,
                "n_frames": 1,
                "totaltime": 0.0,
                "dt": 1.0,
            },
        ),
    ],
)
def test_extract_parameterize(topology, trajectory, reference):
    if trajectory is not None:
        universe = mda.Universe(topology, trajectory)
    else:
        universe = mda.Universe(topology)
    udata = info.extract(universe)

    for key in reference:
        assert udata[key] == pytest.approx(reference[key])

def test_extract_ugmx():
    """Test extract() function for the universe with the standard
    GROMACS test trajectory.
    """
    reference = {
        "n_atoms": 47681,
        "Lx": 80.017006,
        "Ly": 80.017006,
        "Lz": 80.017006,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 90.0,
        "n_frames": 10,
        "totaltime": 900.0000686645508,
        "dt": 100.00000762939453,
    }
    universe = mda.Universe(data.TPR, data.XTC)

    # run the extract() function
    udata = info.extract(universe)

    # assert udata == reference

    # Floating point number comparisons have to be done carefully;
    # https://docs.pytest.org/en/4.6.x/reference.html#pytest-approx
    # works also for integers

    for key in reference:
        assert udata[key] == pytest.approx(reference[key])

### Testing expected failure.
def test_extract_non_universe():
    universe = "This is not a universe."
    with pytest.raises(AttributeError): # this basicaly makes failing a pass
        udata = info.extract(universe)


def test_extract_fixture(u_gmx):
    universe, reference = u_gmx
    # run the extract() function
    udata = info.extract(universe)

    # assert udata == reference

    # Floating point number comparisons have to be done carefully;
    # https://docs.pytest.org/en/4.6.x/reference.html#pytest-approx
    # works also for integers

    for key in reference:
        assert udata[key] == pytest.approx(reference[key])

def test_extract_ucharmm(u_charmm):
    universe, reference = u_charmm
    udata = info.extract(universe)
    for key in reference:
        assert udata[key] == pytest.approx(reference[key])


# test where we capture stdout via capsys
# https://docs.pytest.org/en/4.6.x/capture.html
def test_summary(u_gmx, u_charmm, u_pdb, capsys):
    u1 = u_gmx[0]  # just get the universe, ignore the data
    u2 = u_charmm[0]
    u3 = u_pdb[0]
    labels = ["Adk GROMACS", "ADK DIMS", "Adk PDB"]
    info.summary(u1, u2, u3, labels=labels)
    captured = capsys.readouterr()

    assert (
        "+-------------+---------+-----------+-----------+-----------+-------+------+-------+----------+--------------------+--------------------+"
        in captured.out
    )
    for label in labels:
        assert label in captured.out
    assert "47681" in captured.out
    assert "3341" in captured.out
    assert "Lx/Å" in captured.out

