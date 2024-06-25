'''
This is a docstring kinda thing 
=================================
this file contains functions that do things. 
'''

import numpy as np

from MDAnalysis.lib.distances import calc_bonds, minimize_vectors


def minimage_distance(a, b):
    r"""Calculate distances between atoms |r_a - r_b|.

    Compute the distances between corresponding atoms in two groups with the
    same number of atoms :math:`|\mathbf{r}^{a}_{i} - \mathbf{r}^{b}_{i}|` and
    return a numpy array of the distances of shape `(N,)`.

    Periodic boundary conditions are automatically taken into account.

    Parameters
    ----------
    a: AtomGroup
        First group of N atoms.
    b: AtomGroup
        Second group of N atoms.

    Returns
    -------
    distances : numpy.ndarray

    Raises
    ------
    ValueError
        If `a` and `b` have different sizes or if the box dimensions of the
        two groups differ.


    See Also
    --------
    MDAnalysis.lib.distances.calc_bonds

    """

    _check_groups(a, b)

    d = calc_bonds(a, b, box=a.dimensions)
    return d


def minimage_difference_vector(a, b):
    r"""Calculate the minimum image difference vector r_a - r_b

    Compute the distance vector between corresponding atoms
    :math:`\mathbf{r}^{a}_{i} - \mathbf{r}^{b}_{i}` in two groups with the same
    number of atoms and return a numpy array of the distance vectors of shape
    `(N, 3)`.

    Periodic boundary conditions are automatically taken into account.

    Parameters
    ----------
    a: AtomGroup
        First group of N atoms.
    b: AtomGroup
        Second group of N atoms.

    Returns
    -------
    distance_vectors : numpy.ndarray

    Raises
    ------
    ValueError
        If `a` and `b` have different sizes or if the box dimensions of the
        two groups differ.


    See Also
    --------
    MDAnalysis.lib.distances.minimize_vectors

    """

    # ensure that a.dimensions and b.dimensions are the same
    _check_groups(a, b)

    diff_vec = a.positions - b.positions
    if a.dimensions is not None:
        # need minimum image
        r = minimize_vectors(diff_vec, a.dimensions)
    else:
        # no box information, so the difference vector is sufficient
        r = diff_vec
    return r


def _check_groups(a, b):
    """Helper function that checks compatibility of two atomgroups for
    distance calculations

    Raises :exc:`ValueError` if the groups contain different number of atoms or
    have different associated unit cell dimensions.
    """

    if a.n_atoms != b.n_atoms:
        raise ValueError("AtomGroups a and b contain different number of atoms")

    if a.dimensions is None and b.dimensions is None:
        # it's ok when both AtomGroups have no unitcell
        pass
    elif (a.dimensions is None) != (b.dimensions is None):  # xor
        raise ValueError("One AtomGroup does not have unit cell information")
    elif not np.allclose(a.dimensions, b.dimensions):
        raise ValueError("Unit cell dimensions differ between groups.")