'''
This is Information
==================================
this module contains functions that do something. 
'''

def summary(*universe):
    pass

def extract(u):
    """Extract summary information from a universe.

    The information is returned as a dictionary for keys
    - `n_atoms`: number of atoms
    - `Lx`, `Ly`, `Lz`: length of the unit cell in Å (from first frame of the
      trajectory)
    - `alpha`, `beta`, `gamma`: angles of the unit cell in degrees
    - `n_frames`: number of frames in the trajectory
    - `totaltime`: total simulation time in ps
    - `dt`: time between saved frames in the trajectory in ps


    Parameters
    ----------
    u : Universe
        MDAnalysis Universe.

    Returns
    -------
    data : dict

    """

    try:
        Lx, Ly, Lz, alpha, beta, gamma = u.dimensions
    except TypeError:
        # universe without a regular box
        Lx, Ly, Lz, alpha, beta, gamma = 0, 0, 0, 0, 0, 0
    data = {
        "n_atoms": u.atoms.n_atoms,
        "Lx": Lx,
        "Ly": Ly,
        "Lz": Lz,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "n_frames": u.trajectory.n_frames,
        "totaltime": u.trajectory.totaltime,
        "dt": u.trajectory.dt,
    }
    return data

