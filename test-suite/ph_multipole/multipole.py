"""Modules and functions for extracting mutipolar and dielectric tensors."""
# -*- coding: utf-8 -*-
# Copyright (C) 2021-2025 Changpeng Lin
# All rights reserved.

import argparse
import datetime
import itertools
import math
import numbers
import os
import pickle
import warnings
from collections import Counter, OrderedDict, deque
from copy import deepcopy
from io import StringIO

import numpy as np
import scipy.sparse as sparse
import spglib as spg

try:
    from ase import Atoms as baseAtoms

    __ASE__ = True
except ImportError:
    baseAtoms = object
    __ASE__ = False
    warnings.warn(
        "ASE is not installed. Please install ASE <= 3.22.1 for the case when ibrav!=0."
    )


__version__ = "0.1.0"


__logo__ = r"""
 __  __       _ _   _             _
|  \/  |_   _| | |_(_)_ __   ___ | | ___
| |\/| | | | | | __| | '_ \ / _ \| |/ _ \
| |  | | |_| | | |_| | |_) | (_) | |  __/
|_|  |_|\__,_|_|\__|_| .__/ \___/|_|\___|
                     |_|                 """


def welcome(start_time):
    """Print welcome information."""
    print(f"Started on {start_time}" + __logo__ + __version__)


def goodbye(start_time):
    """Print goodbye information and citation."""
    info1 = "You are using `Multipole` program developed in the work:"
    info2 = "Please consider to cite it if you find the program useful."
    authors = "C. Lin, S. Ponce, F. Macheda, F. Mauri, and N. Marzari, "
    journal = "arXiv:2412.18482 (2024)."

    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time

    print(f"Finshed on {end_time}, total time elapsed: {elapsed_time}")
    print("#" * 84)
    print("# " + info1 + " " * 24 + " #")
    print("# " + authors + journal + " #")
    print("# " + info2 + " " * 22 + " #")
    print("#" * 84)


def timeit(job_string=None):
    """A time decorator for measuring time elapsed for a function to run.

    Parameters
    ----------
    job_string : str, optional
        Description of the job to be measured for time elapsed.

    """

    def decorator(func):
        """Decorator function.

        Parameters
        ----------
        func : callable
            Function to be measured for time elapsed.

        """

        def wrapper(*args, **kwargs):
            start_time = datetime.datetime.now()
            result = func(*args, **kwargs)
            end_time = datetime.datetime.now()
            print(job_string + f" finished, time elapsed: {end_time - start_time}")
            return result

        return wrapper

    return decorator


def normalize(X, axis=1):
    """Scale input vectors individually to unit norm.

    This function plays the same role as the one in scikit-learn:
    sklearn.preprocessing.normalize with norm='l2'.

    Parameters
    ----------
    X : numpy.ndarray
        Input matrix.
    axis : int
        Axis along which to normalize. Default is 1.

    Returns
    -------
    numpy.ndarray
        Normalized matrix.

    """
    if axis == 0:
        X = X.T
    X = X / np.linalg.norm(X, axis=1)[:, None]
    if axis == 0:
        X = X.T
    return X


def kron_product(mat, times=2):
    """A wrapper of multiple-time Kronecker product.

    If the input 'times' is smaller than 5, 'kron_product_einsum' will
    be called to calculate the result; Otherwise, 'kron_product_sparse'
    will be used to compute the result.

    Parameters
    ----------
    mat : numpy.ndarray
        Input matrix.
    times : int
        The number of times of Kronecker product to be performed.

    Returns
    -------
    scipy.sparse.coo_matrix
        The resulting Kronecker product in COO sparse form.

    """
    if times < 5:
        result = kron_product_einsum(mat, times)
        return sparse.coo_matrix(result.reshape(3**times, 3**times))
    else:
        return kron_product_sparse(mat, times)


def kron_product_einsum(mat, times=2):
    """Perform multiple-time Kronecker product using einsum.

    Parameters
    ----------
    mat : numpy.ndarray
        Input matrix.
    times : int
        The number of times of Kronecker product to be performed.

    Returns
    -------
    np.ndarray
        The resulting Kronecker product as a tensor.

    """
    input = []
    for i in range(times):
        input.append(mat)
        input.append([i, times + i])

    return np.einsum(*input)


def kron_product_sparse(mat, times=2):
    """Sparse version of multiple-time Kronecker product.

    Parameters
    ----------
    mat : numpy.ndarray
        Input matrix.
    times : int
        The number of times of Kronecker product to be performed.

    Returns
    -------
    scipy.sparse.coo_matrix
        The resulting Kronecker product in COO sparse form.

    """
    result = sparse.kron(mat, mat)
    for n in range(times - 2):
        result = sparse.kron(result, mat)

    return result.tocoo()


def get_permutation_tensor(index_orginal, index_permuted):
    """Compute permutation tensor according to input atom index.

    This function is similar to 'get_permutation_matrix'. However,
    it uses tensor notation and the resulting permutation for IFCs
    is returned as a multi-dimension tensor by using numpy.einsum.

    Parameters
    ----------
    index_original : list or numpy.ndarray
        Original atom index from the input. If it is not given,
        the original atom index will be taken from object itself.
    index_permuted : list or numpy.ndarray
        Atom index of the cluster after permutation.

    Returns
    -------
    numpy.ndarray
        The permutation tensor with the IFC-order related shape.

    """
    order = len(index_orginal)
    Rmat = kron_product_einsum(np.identity(3), order)
    if index_orginal == index_permuted:
        return Rmat
    else:
        input = list(range(order * 2))
        input_org = list(range(order * 2))
        index_tmp = deepcopy(index_permuted)
        for i, idx in enumerate(index_orginal):
            j = index_tmp.index(idx)
            input[order + j] = input_org[order + i]
            index_tmp[j] = -1
        return np.einsum(Rmat, input)


def get_rotation_cartesian(rotmat, cell, eps=1e-4):
    """Convert rotation matrix from crystal to Cartesian coordinate.

    Parameters
    ----------
    rotmat : numpy.ndarray
        (3,3) Rotation matrix in crystal coordinate.
    cell : numpy.ndarray
        (3,3) Lattice vectors of supercell.

    Returns
    -------
    numpy.ndarray
        (3,3) Rotation matrix in Cartesian coordinate.

    """
    crotmat = np.linalg.multi_dot((cell.T, rotmat, np.linalg.inv(cell.T)))
    crotmat = np.where(abs(crotmat) < eps, 0.0, crotmat)
    return crotmat


def block_diag_sparse(mats, order, format=None, dtype=None):
    """Build a block diagonal sparse matrix from provided matrices.

    This function is originally defined in scipy.sparse.block_diag
    which has been modified here to tackle with the special case of
    making a block diagonal sparse matrix from a list of sub-null-space
    for each representative clusters, i.e. a sub-null-space having zero
    vectors.

    https://github.com/scipy/scipy/blob/v1.8.0/scipy/sparse/_construct.py

    Parameters
    ----------
    mats : sequence of matrices
        Input matrices.
    order : int
        The order of IFCs or cluster.
    format : str, optional
        The sparse format of the result (e.g., "csr"). If not given, the matrix
        is returned in "coo" format.
    dtype : dtype specifier, optional
        The data-type of the output matrix. If not given, the dtype is
        determined from that of `blocks`.

    Returns
    -------
    res : sparse matrix

    """
    row = []
    col = []
    data = []
    r_idx = 0
    c_idx = 0
    fc_row = np.power(3, order)
    for a in mats:
        if isinstance(a, (list, numbers.Number)):
            a = sparse.coo_matrix(a)
        if a.shape == (1, 0) or a.shape == (0,):
            nrows, ncols = fc_row, 0
        else:
            nrows, ncols = a.shape
        if sparse.issparse(a):
            a = a.tocoo()
            row.append(a.row + r_idx)
            col.append(a.col + c_idx)
            data.append(a.data)
        else:
            if a.shape != (1, 0) or a.shape != (0,):
                a_row, a_col = np.divmod(np.arange(nrows * ncols), ncols)
                row.append(a_row + r_idx)
                col.append(a_col + c_idx)
                data.append(a.ravel())
        r_idx += nrows
        c_idx += ncols
    row = np.concatenate(row)
    col = np.concatenate(col)
    data = np.concatenate(data)
    return sparse.coo_matrix(
        (data, (row, col)), shape=(r_idx, c_idx), dtype=dtype
    ).asformat(format)


def rref_dense(mat, eps=1e-4):
    """A dense version of full-pivot Gauss-Jordan elimination.

    Parameters
    ----------
    mat : numpy.ndarray
        Input matrix.
    eps : float
        Numeric tolerance.

    Returns
    -------
    mat : numpy.ndarray
        Input matrix in reduced row echelon form.
    pivots : list
        A list of pivots for the reduced row echelon matrix.

    """
    pivots = np.array([], dtype=int).reshape(0, 2)
    mat = mat[np.where(np.linalg.norm(mat, axis=1) > eps)]
    if mat.size == 0:
        return (mat, pivots)

    while True:
        mat_zero = mat.copy()
        if len(pivots) > 0:
            mat_zero[pivots[:, 0]] = 0
            mat_zero[:, pivots[:, 1]] = 0
        if np.max(np.absolute(mat_zero)) == 0:
            break

        row, col = np.unravel_index(np.absolute(mat_zero).argmax(), mat.shape)
        pivot = mat[(row, col)]
        mat_row = mat[row] / pivot
        mat -= np.outer(mat[:, col], mat[row]) / pivot
        mat[row] = mat_row

        pivots = np.append(pivots, np.array([[row, col]]), axis=0)
        zero_rows = np.where(np.linalg.norm(mat, axis=1) < eps)[0]
        if zero_rows.size != 0:
            mat = np.delete(mat, zero_rows, axis=0)
            move = np.sum(zero_rows[:, None] < pivots[:, 0], axis=0)
            pivots[:, 0] -= move

    return (mat, pivots)


def null_space_dense(mat, eps=1e-4):
    """Calculate the null space of the input matrix.

    Parameters
    ----------
    mat : numpy.ndarray
        Input matrix.
    eps : float, optional
        Numeric tolerance. By default 1E-4.

    Returns
    -------
    numpy.ndarray
        Null space of the input matrix.

    """
    rref, pivots = rref_dense(mat, eps)
    shape = rref.shape
    if len(pivots) == 0:
        null_space = np.eye(shape[1])
    else:
        pivots = np.array(pivots)
        col = np.arange(shape[1], dtype="int")
        col = col[np.isin(col, pivots[:, 1], invert=True)]
        if col.shape == (0,):
            null_space = np.array([])
        else:
            null_space = np.zeros((shape[1], col.shape[0]))
            null_space[col, range(col.shape[0])] = 1.0
            null_space[pivots[:, 1, None], range(col.shape[0])] = -rref[
                pivots[:, 0, None], col
            ]

    return null_space


def wavevector_perturb_simple(nq, q_max, q_min=0.01, q_dir=None, rand_seed=None):
    """Randomly add a phonon perturbation.

    It will add a phonon q along a random direction with the magnitude
    given in between (q_min, q_max).

    Parameters
    ----------
    nq : int
        Number of q-points.
    q_max : float
        Maximum modulus of q-point in reduced coordinate, namely Cartesian
        coordinate divided by 2pi/alat.
    q_min : float, optional
        Minimum modulus of q-point in reduced coordinate. By default 0.01.
    q_dir : array_like, optional
            A (3,) vector in Cartesian coordinate representing the
            direction of q-points, i.e. along a high-symmetry line.
            By default None.
    rand_seed : int, optional
        Seed for random number generator. By default None.

    Returns
    -------
    numpy.ndarray
        The perturbing wavevector in shape (nq, 3).

    """
    if q_dir is None:
        rng = np.random.default_rng(rand_seed)
        q_dirs = rng.uniform(-1, 1, (nq, 3))
    else:
        q_dirs = np.tile(q_dir, (nq, 1))
    q_norm = np.linspace(q_min, q_max, nq)
    q_vecs = (
        q_norm.reshape(-1, 1) * q_dirs / np.linalg.norm(q_dirs, axis=1).reshape(-1, 1)
    )

    return q_vecs


def wavevector_perturb_uniform(grid, q_step=0.01):
    """Generate a uniform grid of q-points.

    This returns a grid around the origin with the step size
    of q_step in crystal coordinate.

    Parameters
    ----------
    grid : array_like
        The number of q-points in each direction.
    q_step : float
        The step size of q-points in crystal coordinate.
        By default 0.01.

    Returns
    -------
    numpy.ndarray
        The q-points in crystal coordinate in shape (nq, 3).

    """
    q_pos = []
    for i in range(3):
        q_min = -q_step * (grid[i] - 1) / 2
        q_max = q_step * (grid[i] - 1) / 2
        q_arr = np.linspace(q_min, q_max, grid[i])
        q_pos.append(q_arr)
    q_vecs = np.array(np.meshgrid(*q_pos)).T.reshape(-1, 3)
    q_vecs = q_vecs[~np.all(q_vecs == 0, axis=1)]

    return q_vecs


def get_irreducible_qpoints(qpts, symops, is_time_reversal=True, eps=1e-5):
    """Get irreducible q-points using symmetry operations.

    Parameters
    ----------
    qpts : numpy.ndarray
        The q-points in crystal coordinate.
    symops : dict
        The symmetry operations in crystal coordinate.
    is_time_reversal : bool, optional
        Whether time reversal symmetry is present. By default True.
    eps : float, optional
        The numerical tolerance. By default 1E-5.

    Returns
    -------
    numpy.ndarray
        The irreducible q-points in crystal coordinate.

    """
    q_ir = []
    inds_eq = set()
    nqts = qpts.shape[0]
    for iq in range(nqts):
        if iq in inds_eq:
            continue
        q = qpts[iq]
        qs_rot = np.dot(symops["rotations"], q)
        if is_time_reversal:
            qs_rot = np.vstack((qs_rot, -qs_rot))
        qs_rot = np.unique(qs_rot, axis=0)
        inds = [
            ind[0] if len(ind) > 0 else -1
            for q in qs_rot
            for ind in [
                np.where(np.linalg.norm(qpts - q - np.rint(qpts - q), axis=1) < eps)[0]
            ]
        ]
        inds_eq = inds_eq.union(set(inds))
        q_ir.append(q)

    return np.array(q_ir)


def build_multipole_wavevector_matrix(multipole_space, natom, q_vec):
    """Construct sensing matrix of given wavevector for multipole expansion.

    Parameters
    ----------
    multipole_space : MultipoleSpace
        A space containing all of symmetry-distinct
        representative onsite clusters for multipole
        tensors followed by their orbits.
    natom : int
        The number of atoms in (primitive) unit cell.
    q_vec : numpy.ndarray
        A q-point in Cartesian coordinate.

    Returns
    -------
    numpy.ndarray
        Sensing matrix of input q wavevector.

    """
    max_order = multipole_space.max_order
    clusters = multipole_space.get_multipole_space()

    sensing_mat_block = deque()
    for order in range(2, max_order + 1):
        shape = [3] * order
        size = np.power(3, order)
        col = np.arange(size, dtype=int)
        for orbit in clusters:
            block = np.zeros((natom * 3, size), dtype=complex)
            for cluster in orbit[1:]:
                block_tmp = np.zeros((natom * 3, size), dtype=complex)
                comp_list = np.array(list(np.ndindex(*shape)))
                if order == 2:
                    row = (comp_list + cluster.atom_index * 3)[:, 1]
                    comp_list = np.delete(comp_list, 1, 1)
                else:
                    row = (comp_list + cluster.atom_index * 3)[:, 0]
                    comp_list = np.delete(comp_list, 0, 1)
                q_prods = np.prod(np.take(q_vec, comp_list), axis=1)
                np.add.at(block_tmp, (row, col), q_prods)
                Gamma = cluster.get_crotation_tensor(order).toarray()
                block += (
                    block_tmp.dot(Gamma)
                    / math.factorial(order - 1)
                    * (-1j) ** (order - 1)
                )
            sensing_mat_block.append(block)

    return np.hstack(sensing_mat_block)


def build_dielectric_wavevector_matrix(epsil_order, q_vec):
    """Construct sensing matrix for expanding dielectric function.

    Parameters
    ----------
    epsil_order : int
        The maximum order of dielectric tensor used to expand
        macroscopic dielectric function.
    q_vec : numpy.ndarray
        A q-point in Cartesian coordinate.

    Returns
    -------
    numpy.ndarray
        Sensing matrix of input q wavevector.

    """
    sensing_mat_block = deque()
    for order in range(2, epsil_order + 1, 2):
        shape = [3] * order
        comp_list = np.array(list(np.ndindex(*shape)))
        block = np.prod(np.take(q_vec, comp_list), axis=1)
        sensing_mat_block.append(block)

    return np.hstack(sensing_mat_block)


def read_charge_density_response(natom, filename="drho.dat"):
    """Read charge density response due to a phonon perturbation.

    This function can be also used to read effective charges.

    Parameters
    ----------
    natom : int
        The number of atoms in (primitive) unit cell.
    filename : str
        The filename of charge density response.

    Returns
    -------
    numpy.ndarray
        Charge density response of G=0 component in the shape
        (natom, 3).

    """
    # drho_tmp = np.loadtxt(filename)[:natom, 1:].reshape([natom, 3, 2])
    drho_tmp = np.loadtxt(filename)[:natom, :].reshape([natom, 3, 2])
    drho = np.empty([natom, 3], dtype=complex)
    drho.real = drho_tmp[:, :, 0]
    drho.imag = drho_tmp[:, :, 1]

    return drho


def read_inverse_dielectric_function(filename="drhodv.dat"):
    """Read inverse dielectric function at a q-point.

    Parameters
    ----------
    filename : str
        The filename of drhodv calculated by ph.x of
        Quantum ESPRESSO.

    Returns
    -------
    numpy.ndarray
        A vector containing epsm1, 1/epsm1, chi and 1+vq*chi
        in sequence.

    """
    data = np.loadtxt(filename, usecols=(0, 1))
    drhodv = np.empty(data.shape[0], dtype=complex)
    drhodv.real = data[:, 0]
    drhodv.imag = data[:, 1]

    return drhodv


def write_ph_input(q_vec, filename, ph_header="ph.in"):
    """Write ph.x input file with the given q-point.

    Parameters
    ----------
    q_vec : numpy.ndarray (3,)
        A q-point in reduced coordinate.
    filename : str
        The filename for ph.x input to be written.
    ph_header : str
        The filename of ph.x input header.

    """

    iq = filename.split(".")[-1]
    with open(ph_header) as fd:
        header = fd.readlines()
    ph_in = StringIO()
    for line in header:
        if "{}" in line:
            ph_in.write(line.format(iq))
        else:
            ph_in.write(line)
    ph_in.write("{0[0]:>13.8f} {0[1]:>13.8f} {0[2]:>13.8f}\n".format(q_vec.tolist()))
    with open(filename, "w") as finalf:
        finalf.write(ph_in.getvalue())
    ph_in.close()


def write_symops(symops, filename="spglib.symops"):
    """Write symmetry operations searched by spglib into file."""
    nsym = symops["translations"].shape[0]
    with open(filename, "w") as f:
        f.write("%d\n" % nsym)
        for i in range(nsym):
            f.write(
                "{0[0][0]:4d} {0[0][1]:4d} {0[0][2]:4d} {0[1][0]:4d} {0[1][1]:4d} {0[1][2]:4d} {0[2][0]:4d} {0[2][1]:4d} {0[2][2]:4d}\n".format(
                    symops["rotations"][i].tolist()
                )
            )
            f.write(
                "{0[0]:8.6f} {0[1]:8.6f} {0[2]:8.6f}\n".format(
                    symops["translations"][i].tolist()
                )
            )


class InputParser(argparse.ArgumentParser):
    """Multipole input parser inheriting from argparse.ArgumentParser class.

    Parameter settings from user-defined command line for manipulating
    Multipole program is realized in this class. Input arguments that user
    can set must appear in upper case, while flags and commands must
    appear in lower case.

    https://docs.python.org/3/library/argparse.html

    """

    def __init__(self):
        """Initialize default settings."""
        super(InputParser, self).__init__()
        self.prog = "MULTIPOLE " + __version__
        self.description = "A program for extracting multipole and dielectric tensors."
        self.usage = __logo__ + __version__ + "\n" + "multipole.py -x ... [--flags ...]"
        self.epilog = "MULTIPOLE program command-line inteface. Enjoy !^_^!"
        self.formatter_class = argparse.ArgumentDefaultsHelpFormatter

        """Get default settings"""
        self.init_args()
        self.defaults = self.parse_args([])

        """Get user-defined settings"""
        settings = deepcopy(self.defaults)
        self.settings = self.parse_args(namespace=settings)

    def init_args(self):
        """Initialize arguments and help message."""

        self.add_argument(
            "-c",
            "--cell",
            dest="CELL_FILENAME",
            action="store",
            default="scf.in",
            type=str,
            help="""Primitive cell filename.""",
        )
        self.add_argument(
            "--alat",
            dest="ALAT",
            action="store",
            default=None,
            type=str,
            help="""The celldm(1) or A parameter in pwscf input, required when ASE
                    is present. The format must be like `1.0A` or `1.0B` where `A`
                    and `B` are the units angstrom and bohr, respectively.""",
        )
        self.add_argument(
            "--cry_basis",
            dest="CRY_BASIS",
            action="store_true",
            default=False,
            help="""True to deal with rotation matrix in crystal coordinate (not
                    implemented yet). Please keep it False in Cartesian coordinate.""",
        )
        self.add_argument(
            "-e",
            "--mpe",
            dest="MP_EXPAND",
            action="store_true",
            default=False,
            help="""This option asks code to perform multipole expansion.
                    If `epsil_order` is set, dielectric expansion will be
                    also performed.""",
        )
        self.add_argument(
            "--eps",
            dest="EPS",
            action="store",
            default=1e-4,
            type=float,
            help="Numerical tolerance of matrix norm for constructing null space.",
        )
        self.add_argument(
            "--epsil_kernel",
            dest="EPSIL_KERNEL",
            action="store",
            default=1,
            type=int,
            choices=[1, 2],
            help=r"""The kernel for dielectric screening function.
                    1 for $1+v(q)*\chi(q)$ and 2 for $1/epsilon(q)^{-1}$.
                    This option is only available when the macroscopic potential is
                    not removed, i.e. 'lmacro=.FALSE.' in DFPT, otherwise the kernel
                    $1/epsilon(q)^{-1}$ is used.""",
        )
        self.add_argument(
            "--epsil_order",
            dest="EPSIL_ORDER",
            action="store",
            default=None,
            type=int,
            help="""The highest order used to perform the dielectric expansion.
                    This value must be an even integer.""",
        )
        self.add_argument(
            "-f",
            dest="FIT",
            action="store_true",
            default=False,
            help="This option asks code to fit multipole and dielectric expansions.",
        )
        self.add_argument(
            "--fix_order",
            dest="FIX_ORDER",
            action="store",
            default=None,
            type=int,
            help="""If an order is set, the code will read multipole
                    tensors up to the specified order from files and
                    fix them during fitting procedure.""",
        )
        self.add_argument(
            "--fix_epsil_order",
            dest="FIX_EPSIL_ORDER",
            action="store",
            default=None,
            type=int,
            help="""If an order is set, the code will read dielectric
                    tensors up to the specified order from files and
                    fix them during fitting procedure.""",
        )
        self.add_argument(
            "--ir_q",
            dest="IR_Q",
            action="store_true",
            default=False,
            help="Set to True to use only irreducible q-points for ph.x.",
        )
        self.add_argument(
            "--mesh",
            dest="MESH",
            action="store",
            nargs=3,
            default=None,
            type=int,
            help="Uniform mesh grid for multipole and dielectric expansions.",
        )
        self.add_argument(
            "--mesh_step",
            dest="MESH_STEP",
            action="store",
            default=0.01,
            type=float,
            help="""A wavenumber step in crystal coordinate, only used
                    together with `mesh`.""",
        )
        self.add_argument(
            "--order",
            dest="MP_ORDER",
            action="store",
            default=2,
            type=int,
            help="The Highest order included in multipole expansion.",
        )
        self.add_argument(
            "-n",
            "--nq",
            dest="NQPTS",
            action="store",
            default=None,
            type=int,
            help="Number of q-points for charge density response calculations.",
        )
        self.add_argument(
            "--non_polar",
            dest="NON_POLAR",
            action="store_true",
            default=False,
            help="""Set Born effective charge tensors to zero in the case of
                    non-polar solids.""",
        )
        self.add_argument(
            "-p",
            "--sensing",
            dest="SENSING_MAT",
            action="store_true",
            default=False,
            help="""With this option, the code will generate q-points for
                    charge density response calculation in ph.x and construct
                    the sensing matrix for multipole (and dielectric) expansion.""",
        )
        self.add_argument(
            "--pred",
            dest="PREDICT",
            action="store_true",
            default=False,
            help="""This option asks code to predict charge density response,
                    bare effective charges and dielectric function from
                    multipole and dielectric expansions.""",
        )
        self.add_argument(
            "--q_dir",
            dest="Q_DIR",
            action="store",
            nargs=3,
            default=None,
            type=int,
            help="""The q-point line direction in crystal coordinate. If not set,
                    the direcion for each q-point will be generated randomly.""",
        )
        self.add_argument(
            "--q_file",
            dest="Q_FILE",
            action="store_true",
            default=False,
            help="""If True, read q-points from file.""",
        )
        self.add_argument(
            "--q_min",
            dest="Q_MIN",
            action="store",
            default=0.01,
            type=float,
            help="Minimum magnitude of wavenumber in reduced coordinate.",
        )
        self.add_argument(
            "--q_max",
            dest="Q_MAX",
            action="store",
            default=0.05,
            type=float,
            help="Maximum magnitude of wavenumber in reduced coordinate.",
        )
        self.add_argument(
            "--read",
            dest="MP_READ",
            action="store_true",
            default=False,
            help="""If True, read multipole and dielectric tensors from files when
                    predicting charge density response, bare effective charges and
                    dielectric function using multipole expansion.""",
        )
        self.add_argument(
            "--screened",
            dest="SCREENED",
            action="store_true",
            default=False,
            help="""Set to True if the charge density response is in the screened
                    type, i.e. the macroscopic potential is not removed in DFPT
                    calculations.""",
        )
        self.add_argument(
            "--seed",
            dest="RAND_SEED",
            action="store",
            default=None,
            type=int,
            help="""Seed for pseudo random number generator. This is
                    useful for reproducing the same results, e.g. fixing
                    the randomly generated perturbations.""",
        )
        self.add_argument(
            "--symprec",
            dest="SYMPREC",
            action="store",
            default=1e-5,
            type=float,
            help="Tolerance for finding symmetry using spglib.",
        )
        self.add_argument(
            "-v",
            "--version",
            action="version",
            version="%(prog)s",
            help="Show version information.",
        )
        self.add_argument(
            "--write_symops",
            dest="WRITE_SYMOPS",
            action="store_true",
            default=False,
            help="""True to write space group symmetry operations into file.""",
        )


class Atoms(baseAtoms):
    """Adapted version of Atoms class inheriting from aseAtoms class.

    https://gitlab.com/ase/ase/-/blob/master/ase/atoms.py

    """

    """Constants"""
    ANGSTROM_TO_BOHR = 1.8897261246257702

    def __init__(
        self,
        aseAtoms=None,
        alat=None,
        cell=None,
        positions=None,
        scaled_positions=None,
        symbols=None,
        numbers=None,
    ):
        """Initialize Atoms object."

        Parameters
        ----------
        aseAtoms : ase.Atoms
            An ASE Atoms object. Default is None.
        alat : float
            Lattice constant in pwscf. Default is None.
        cell : numpy.ndarray
            Lattice vectors of primitive cell. Default is None.
        positions : numpy.ndarray
            Positions of atoms in Cartesian coordinate in the shape of (natom, 3).
            Default is None.
        scaled_positions : numpy.ndarray
            Positions of atoms in crystal coordinate in the shape of (natom, 3).
            Default is None.
        symbols : list
            Chemical symbols of atoms in the shape of (natom,). Default is None.
        numbers : list
            Atom type index of atoms in the shape of (natom,). Default is None.

        """
        if aseAtoms is not None:
            super(Atoms, self).__init__(aseAtoms)
        else:
            self._cell = cell
            self._positions = positions
            self._scaled_positions = scaled_positions
            self._symbols = symbols
            self._numbers = numbers
        self._alat = alat

    @property
    def cell(self):
        """numpy.ndarray : lattice vectors of primitive cell."""
        if __ASE__:
            return super().cell.real
        else:
            return self._cell.copy()

    @property
    def positions(self):
        """numpy.ndarray : positions of atoms in Cartesian coordinate."""
        if __ASE__:
            return super().positions
        else:
            return self._positions.copy()

    @property
    def scaled_positions(self):
        """numpy.ndarray : positions of atoms in crystal coordinate."""
        if __ASE__:
            return super().positions.dot(np.linalg.inv(self.cell))
        else:
            return self._scaled_positions.copy()

    def reciprocal(self):
        """numpy.ndarray : reciprocal lattice vectors."""
        if __ASE__:
            return super().cell.reciprocal().real
        else:
            return np.linalg.inv(self.cell).T

    def get_global_number_of_atoms(self):
        """Return the global number of atoms."""
        if __ASE__:
            return super().get_global_number_of_atoms()
        else:
            return len(self._numbers)

    def get_chemical_formula(self):
        """Return the chemical formula."""
        if __ASE__:
            return super().get_chemical_formula()
        else:
            formula = ""
            for symbol, count in Counter(self._symbols).items():
                if count == 1:
                    formula += symbol
                else:
                    formula += symbol + str(count)
            return formula

    def get_chemical_symbols(self):
        """Return the chemical symbols."""
        if __ASE__:
            return super().get_chemical_symbols()
        else:
            return self._symbols.copy()

    def get_atomic_numbers(self):
        """Return the atomic numbers."""
        if __ASE__:
            return super().get_atomic_numbers()
        else:
            return self._numbers.copy()

    @property
    def symops(self):
        """dict : symmetry operations from space group."""
        return self._symops

    @classmethod
    def read(cls, filename="scf.in", format="espresso-in", alat=None):
        """Read structure from file.

        Note that if ASE is not installed, the structure will be read
        using the built-in `read` method of Atoms class which only
        supports reading from Quantum ESPRESSO input file with ibrav=0.

        Parameters
        ----------
        filename : str
            Name of file containing structure. Default is "scf.in".
        format : str
            Format of file. Default is "espresso-in".
        alat : str
            Lattice constant in pwscf. Default is None.

        Returns
        -------
        Atoms
            Structure object.

        """
        if __ASE__:
            import ase.io as aseio

            if alat is None:
                warnings.warn(
                    "Lattice parameter `alat` is not set, q-point unit maybe wrong."
                )
            else:
                if "B" in alat.upper():
                    alat = float(alat.upper().split("B")[0]) / cls.ANGSTROM_TO_BOHR
                elif "A" in alat.upper():
                    alat = float(alat.upper().split("A")[0])
                else:
                    raise ValueError("Unknown lattice constant unit.")
            cell = aseio.read(filename, format=format)
            return cls(aseAtoms=cell, alat=alat)
        else:
            return cls.read_pwscf(filename)

    def get_space_group(self, symprec=1e-5, angle_tolerance=-1.0, symbol_type=0):
        """Return space group information.

        Parameters
        ----------
        symprec : float
            Tolerance for determining symmetry. Default is 1e-5.
        angle_tolerance : float
            Tolerance for determing angle of cell. Default is -1.0.
        symbol_type : int
            Space group symbol type: 0 for international symbol and
            1 for schoenflies symbol. Default is 0.

        Returns
        -------
        str
            Space group symbol and number as a string.

        """
        return spg.get_spacegroup(
            self.to_spglib_tuple(),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            symbol_type=symbol_type,
        )

    def get_symmetry(self, symprec=1e-5, angle_tolerance=-1.0):
        """Find and return symmetry operations.

        Parameters
        ----------
        symprec : float
            Tolerance for determining symmetry. Default is 1e-5.
        angle_tolerance : float
            Tolerance for determing angle of cell. Default is -1.0.

        Returns
        -------
        dict
            A dictionary containing symmetry operations as key-value pairs.
            The keys are as follows:
            "rotations": rotation matrices
            "translations": translation vectors
            "equivalent_atoms": equivalent atoms

        """
        if not hasattr(self, "_symops"):
            self._symops = spg.get_symmetry(
                self.to_spglib_tuple(),
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
            self._nsym = self._symops["translations"].shape[0]
        return self._symops

    def get_number_of_symmetries(self):
        """Return number of space group symmetry operations."""
        if not hasattr(self, "_nsym"):
            self.get_symmetry()
        return self._nsym

    def get_lattice_constant(self, unit="Angstrom"):
        """Calculate lattice constant of the first lattice vector.

        Parameters
        ----------
        unit : str
            The unit of lattice constant to be computed,
            can be "Angstrom" or "Bohr".

        Returns
        -------
        float
            Lattice constant of the first lattice vector.

        """
        if hasattr(self, "_alat") and self._alat is not None:
            alat = self._alat
        else:
            alat = np.linalg.norm(self.cell[0])
        if unit.upper() == "ANGSTROM":
            return alat
        elif unit.upper() == "BOHR":
            return alat * self.ANGSTROM_TO_BOHR
        else:
            raise ValueError(f"Unknown unit: {unit}")

    def get_volume(self, unit="Angstrom"):
        """Calculate the volume of cell.

        Parameters
        ----------
        unit : str
            The unit of lattice constants to be computed,
            can be "Angstrom" or "Bohr".

        Returns
        -------
        Float
            Volume of cell.

        """
        volume = abs(self.cell[0].dot(np.cross(self.cell[1], self.cell[2])))
        if unit.upper() == "ANGSTROM":
            return volume
        elif unit.upper() == "BOHR":
            return volume * self.ANGSTROM_TO_BOHR**3
        else:
            raise ValueError(f"Unknown unit: {unit}")

    def to_spglib_tuple(self):
        """Convert and return structure information as a tuple.

        Returns
        -------
        tuple
            The tuple consists of
            (lattice vectors,
             scaled positions of atoms in the cell,
             the corresponding atom types).

        """
        return (self.cell, self.scaled_positions, self.get_atomic_numbers())

    @classmethod
    def read_pwscf(cls, filename):
        """Read structure from Quantum ESPRESSO input file.

        Parameters
        ----------
        filename : str
            Name of file containing structure.

        Returns
        -------
        Atoms
            Crystal structure object.

        """
        with open(filename, "r") as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines]
            lines = [line.replace(",", "") for line in lines]
            lines = [line.replace("(", "") for line in lines]
            lines = [line.replace(")", "") for line in lines]
            lines = list(filter(None, lines))

        is_pwscf = False
        for n, line in enumerate(lines):
            line_split = line.split()
            if line_split[0].lower() == "&system":
                is_pwscf = True
            if "=" in line_split:
                if line_split[0].lower() == "ibrav":
                    ibrav = int(line_split[2])
                    if ibrav != 0:
                        raise ValueError(
                            "Only ibrav=0 is supported if ASE is not installed."
                        )
                if line_split[0].lower() == "nat":
                    nat = int(line_split[2])
                if line_split[0].lower() == "ntyp":
                    ntyp = int(line_split[2])
                if line_split[0].lower() == "celldm1":
                    alat = float(line_split[2]) / cls.ANGSTROM_TO_BOHR
                if line_split[0].upper() == "A":
                    alat = float(line_split[2])
            if line_split[0].upper() == "CELL_PARAMETERS":
                cell = np.zeros((3, 3))
                for i in range(3):
                    cell[i] = list(map(float, lines[n + i + 1].split()))
                if "bohr" in line_split:
                    cell /= cls.ANGSTROM_TO_BOHR
                    alat = np.linalg.norm(cell[0])
                elif "angstrom" in line_split:
                    alat = np.linalg.norm(cell[0])
                else:
                    cell *= alat
            if line_split[0].upper() == "ATOMIC_POSITIONS":
                try:
                    cell
                except NameError or UnboundLocalError:
                    raise ValueError(
                        "CELL_PARAMETERS must be set before ATOMIC_POSITIONS."
                    )
                atom_pos = np.zeros((nat, 3))
                symbols = []
                for i in range(nat):
                    atom_info = lines[n + i + 1].split()
                    symbols.append(atom_info[0])
                    atom_pos[i] = list(map(float, atom_info[1:]))
                if "alat" in line_split:
                    positions = atom_pos * alat
                    scaled_positions = positions.dot(np.linalg.inv(cell))
                elif "crystal" in line_split:
                    positions = atom_pos.dot(cell)
                    scaled_positions = atom_pos
                elif "bohr" in line_split:
                    positions = atom_pos / cls.ANGSTROM_TO_BOHR
                    scaled_positions = positions.dot(np.linalg.inv(cell))
                elif "angstrom" in line_split:
                    positions = atom_pos
                    scaled_positions = positions.dot(np.linalg.inv(cell))
                else:
                    raise ValueError("Unknown unit for atomic positions.")
            if line_split[0].upper() == "ATOMIC_SPECIES":
                types = []
                for i in range(ntyp):
                    atom_info = lines[n + i + 1].split()
                    types.append(atom_info[0])
        numbers = []
        for i in range(nat):
            for j in range(ntyp):
                if symbols[i] == types[j]:
                    numbers.append(j)
                    break
        if not is_pwscf:
            raise ValueError("Not a Quantum ESPRESSO input file.")

        return cls(
            alat=alat,
            cell=cell,
            positions=positions,
            scaled_positions=scaled_positions,
            symbols=symbols,
            numbers=numbers,
        )


class Site(object):
    """Class for representing an atom site in crystals.

    It can be also used to represent a onsite cluster used in
    constructing multipole tensors in multipole expansion.

    """

    def __init__(
        self,
        atom_pos,
        atom_index,
        atom_type,
        atom_symbol,
        lrep=True,
        rotmat=None,
        trans=None,
    ):
        """Initialization function.

        Parameters
        ----------
        atom_pos : array_like
            Atom position in crystal coordinate in the shape of (3,).
        atom_index : int
            The corresponding atom index within the unit cell or supercell.
        atom_type : int
            The corresponding atom type or atomic number
        atom_symbol : str
            The corresponding chemical symbol
        lrep : bool
            True for representative site and False for those in orbit.
        rotmat : (3,3) numpy.ndarray
            Rotation matrix in crystal coordinate.
        trans : (3,0) numpy.ndarray
            Translational vector in crystal coordinate.

        """
        self._lrep = lrep
        self._atom_pos = atom_pos
        self._atom_index = atom_index
        self._atom_type = atom_type
        self._atom_symbol = atom_symbol
        if not self._lrep:
            self._rotmat = rotmat
            self._trans = trans

    @property
    def is_rep(self):
        """Check if the site is representative."""
        return self._lrep

    @property
    def atom_pos(self):
        """atom_pos: atomic position for the site."""
        return self._atom_pos

    @property
    def atom_index(self):
        """atom_index: atom index in unit cell or supercell for the site."""
        return self._atom_index

    @property
    def atom_type(self):
        """atom_type: atom type for the site."""
        return self._atom_type

    @property
    def atom_symbol(self):
        """atom_symbol: chemical symbol for the site."""
        return self._atom_symbol

    @property
    def rotmat(self):
        """rotmat: Rotation matrix in crystal coordinate."""
        if self._lrep:
            raise AttributeError(
                "Rotation matrix is not available for representative site."
            )
        return self._rotmat

    @property
    def crotmat(self):
        """crotmat: Rotation matrix in Cartesian coordinate."""
        if self._lrep:
            raise AttributeError(
                "Rotation matrix is not available for representative site."
            )
        return self._crotmat

    @crotmat.setter
    def crotmat(self, crotmat):
        """Set (3,3) rotation matrix in Cartesian coordinate.

        Parameters
        ----------
        crotmat : numpy.ndarray
            (3,3) Rotation matrix in Cartesian coordinate.

        """
        self._crotmat = crotmat

    @property
    def trans(self):
        """trans: Translational vector in crystal coordinate."""
        if self._lrep:
            raise AttributeError(
                "Translational vector is not available for representative site."
            )
        return self._trans

    def get_crotation_tensor(self, order):
        """Return rotation tensor for multipole tensor (Cartesian coordinate).

        Parameters
        ----------
        order : int
            The order of rotation tensor to be calculated

        """
        crot_tensor = kron_product(self._crotmat, order)
        return crot_tensor

    def set_isotropy_group(self, isotropy):
        """Set the isotropy group of a representative site.

        Parameters
        ----------
        isotropy : list(dict)
            Isotropy group of a representative atom site. It is a list
            and each element therein is a dictionary corresponding to
            a cluster that is isotropic to the representative one. The
            dictionary should have keys of 'rotmat', 'crotmat' and 'trans'.

        """
        self._isotropy = isotropy

    def get_isotropy_group(self):
        """Return the isotropy group of a representative cluster."""
        return self._isotropy

    def build_isotropy_symmetry_constraints(self, order, cry_basis=False):
        """Construct symmetry constraints from isotropy group.

        Parameters
        ----------
        order : int
            The order of onsite cluster.
        cry_basis : bool
            True to deal with rotation matrix in crystal coordinate,
            False for Cartesian coordinate.

        Returns
        -------
        list(scipy.sparse.coo_matrix)
            A list of isotropy symmerty constraints in COO sparse matrix form.

        """
        if not self._lrep:
            raise NotImplementedError(
                "Isotropy symmetry cannot be built for a non-representative site or cluster."
            )

        if cry_basis:
            rot_fmt = "rotmat"
        else:
            rot_fmt = "crotmat"
        constraints = []
        size = np.power(3, order)
        for cluster in self._isotropy:
            Gamma = kron_product(cluster[rot_fmt], order)
            constraints.append(Gamma - sparse.identity(size))

        return constraints

    def build_permutation_symmetry_constraints(self, order):
        """Construct permutation symmetry constraints.

        This is only applicable for the onsite cluster which has
        the order larger than 1.

        Parameters
        ----------
        order : int
            The order of Onsite cluster.

        Returns
        -------
        list(scipy.sparse.coo_matrix)
            A list of permutation symmerty constraints in COO sparse matrix form.

        """
        if not self._lrep:
            raise NotImplementedError(
                "Symmetry constraints cannot be built for non-representative cluster."
            )
        if order < 2:
            raise ValueError(
                "Permutation symmetry cannot be built for "
                + f"cluster of order lower than {order}."
            )

        permuted_list = []
        for i in range(1, order - 1):
            for j in range(i + 1, order):
                idx = list(range(order))
                idx[i], idx[j] = idx[j], idx[i]
                permuted_list.append(idx)

        constraints = []
        size = np.power(3, order)
        idx = list(range(order))
        for item in permuted_list:
            Rmat = get_permutation_tensor(idx, item).reshape((size, size))
            cons_mat = sparse.identity(size) - sparse.coo_matrix(Rmat)
            constraints.append(cons_mat)

        return constraints

    def __repr__(self):
        """Return the atom site information."""
        txt = "Site(index: {!r}, atom: {!r}, atom_type: {!r}, is_rep: {!r})"
        return txt.format(self._atom_index, self._atom_pos, self._atom_type, self._lrep)

    def __eq__(self, clus):
        """Check if two atom sites are equal based on atomic index."""
        return sorted(self._atom_index) == sorted(clus.atom_index)

    def __ne__(self, clus):
        """Check if two atom sites are not equal based on atomic index."""
        return sorted(self._atom_index) != sorted(clus.atom_index)


class Optimizer(object):
    """Multipole and dielectric tensor optimizer.

    This class adopts the ordinary least-square to extract mulipole
    and dielectric tensors from a perturbation expansion. It mainly
    implements the `fit` and `predict` methods as well as the method
    to score the model `score_metrics` and `set_model_metrics`.

    """

    def __init__(self):
        """Initialize an optimizer."""
        self.results = {}
        self.metrics = {}

    def fit(self, X, y):
        """Fit regression model with the specified method and parameters.

        Parameters
        ----------
        X : numpy.ndarray
            2D sensing matrix.
        y : numpy.ndarray
            1D interatomic force array.
        weights : numpy.ndarray
            Weights for each sample in F.

        """
        self.coef_ = np.linalg.lstsq(X, y, rcond=None)[0]
        self.n_features_in_ = X.shape[1]

        """Compute various metrics"""
        self.set_model_metrics(X, y)

    def predict(self, X):
        """Predict using the linear model.

        Parameters
        ----------
        X : numpy.ndarray of shape (n_samples, n_features)
            Samples.

        Returns
        -------
        numpy.ndarray
            Predicted values.

        """
        return X.dot(self.coef_)

    def set_model_metrics(self, X, y):
        """Compute various metrics for the model.

        Parameters
        ----------
        X : numpy.ndarray
            2D sensing matrix.
        y : numpy.ndarray
            1D interatomic force array.

        """
        self.results["parameters"] = self.coef_
        self.results["n_featrues"] = self.n_features_in_
        self.results["n_nonzero_featrues"] = np.count_nonzero(self.coef_)

        """Compute metrics"""
        eps = np.finfo(np.float64).eps
        y_pred = self.predict(X)
        y_err = np.abs(y_pred - y)
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)

        # Rank of coefficient matrix
        self.metrics["rank"] = np.linalg.matrix_rank(X)
        # Relative norm error
        self.metrics["rne"] = np.sqrt(np.dot(y_err, y_err) / np.dot(y, y))
        # R^2 score
        self.metrics["r2_score"] = 1 - ss_res / ss_tot
        # Mean squared error
        self.metrics["mse"] = np.mean((y - y_pred) ** 2)
        # Root mean squared error
        self.metrics["rmse"] = np.sqrt(self.metrics["mse"])
        # Mean absolute percentage error
        self.metrics["mape"] = np.mean(np.abs(y_pred - y) / np.maximum(np.abs(y), eps))


class MultipoleSpace(object):
    """Class defining a container-like object for multipole expansion.

    This space is a container that includes all of symmetry-distinct
    atom sites and their orbits (obtained by symmetry operations of
    a space group), which are further represented by class Site.

    """

    """Global configuration for multipole expansion."""
    MultipoleSpaceFile = "ms.pkl"

    def __init__(self, cell, symops, mpolar_space, max_order):
        """Called only by class method do_multipole_expansion.

        Parameters
        ----------
        cell : numpy.ndarray
            Lattice vectors of unit cell.
        symops : dict
            Symmetry operations of space group.
            It has the keys "translations" and "rotations".
        mpolar_space : list
            Each element in the list starts by a representative
            site followed by its orbit.
        max_order : int
            highest order in generating cluster space.

        """
        self._cell = cell
        self._symops = symops
        self._MS = mpolar_space
        self._max_order = max_order

    @property
    def max_order(self):
        """max_order: the maximum order of multipole expansion."""
        return self._max_order

    @property
    def cell(self):
        """cell: lattice vectors of unit cell."""
        return self._cell

    @property
    def size(self):
        """Return the number of representative sites in this space."""
        return len(self._MS)

    def get_symmetry(self):
        """Return space group symmetry for (primitive) unit cell."""
        return self._symops

    def get_multipole_space(self):
        """Return the multipole space."""
        return self._MS

    def get_representative_sites(self):
        """Return the representative sites in this space."""
        return [orbit[0].atom_index for orbit in self._MS]

    def get_number_of_components_all_order(self):
        """Return number of multipole tensor components for each order.

        Returns
        -------
        List
            A list of number of multipole tensor components before
            symmetrization for each order.

        """
        ncomp = []
        for order in range(2, self.max_order + 1):
            ncomp.append(np.power(3, order) * self.size)
        return ncomp

    def predit(self, sensing_mat_q1, sensing_mat_q2=None):
        """Predit dielectric properties using mulipole expansion.

        The bare effective charge response to a phonon perburbation
        is calculated based on multipole expansion. If ``sensing_mat_q2``
        is not None, macroscopic dielectric function will be also
        calculated.

        Parameters
        ----------
        sensing_mat_q1 : numpy.ndarray
            2D sensing matrix for charge response.
        sensing_mat_q2 : numpy.ndarray
            2D sensing matrix for dielectric function.

        Returns
        -------
        numpy.ndarray
            Predicted charge response and dielectric function.

        """
        zeff_pred = np.dot(
            sensing_mat_q1.toarray(), self.get_multipole_tensors(lrep=True)
        )
        if sensing_mat_q2 is not None:
            epsil_pred = np.dot(sensing_mat_q2.toarray(), self.get_dielectric_tensors())
            return zeff_pred, epsil_pred
        else:
            return zeff_pred

    def write(self, filename="ms.pkl"):
        """Write MultipoleSpace instance into pickle file.

        # TODO: only dump necessary attributes.

        Parameters
        ----------
        filename : str
            Filename to save multipole space.

        """
        with open(filename, "wb") as fd:
            pickle.dump(self, fd)

    @staticmethod
    def read(filename="ms.pkl"):
        """Read and create MultipoleSpace instance from pickle file.

        Parameters
        ----------
        filename : str
            Filename for MultipoleSpace instance to read from.

        Returns
        -------
        MultipoleSpace object

        """
        with open(filename, "rb") as fd:
            return pickle.load(fd)

    def print_multipole_space_info(self):
        """Print essential info for this multipole space."""
        print("Atoms to be calculated:")
        for item in self._MS:
            idx = item[0].atom_index
            symbol = item[0].atom_symbol
            print(f"- Index {idx} | {symbol}")

    @classmethod
    def run(cls, settings, pcell):
        """Run multipole expansion for charge response and dielectric function.

        This is a class method that creates a MultipoleSpace instance and should
        be only called by the main program.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.
        pcell : Pheasy.Atoms
            Primitive unit cell.

        Returns
        -------
        MultipoleSpace
            An instance of class MultipoleSpace containing all the symmetry-
            distinct sites for multipoles followed by their orbits.

        """
        max_order = settings.MP_ORDER
        epsil_order = settings.EPSIL_ORDER
        MP_space = None

        if max_order == 2:
            log_info = "dynamical Born effective charges."
        elif max_order == 3:
            log_info = "dynamical quadrupoles."
        elif max_order == 4:
            log_info = "dynamical octupoles."
        elif max_order == 5:
            log_info = "dynamical hexadecapoles."
        elif max_order == 6:
            log_info = "dynamical triacontadipoles."
        else:
            log_info = f"{max_order}-order."

        if settings.MP_EXPAND:
            print("Performing multipole expansion up to " + log_info)
            if epsil_order is not None:
                print(
                    f"The {epsil_order}-order expansion "
                    + "of dielectric function is included."
                )

            MP_space = cls.do_multipole_expansion(pcell, max_order)
            MP_space.write(cls.MultipoleSpaceFile)
        else:
            """Read multipole space from file."""
            if os.path.isfile(cls.MultipoleSpaceFile):
                MP_space = cls.read(cls.MultipoleSpaceFile)
                print(
                    "Reading and generating multipole space from file up to " + log_info
                )
                MP_space.print_multipole_space_info()

        return MP_space

    @staticmethod
    @timeit("Multipole expansion")
    def do_multipole_expansion(pcell, max_order, eps=1e-4):
        """Find symmetry-distinct sites and their orbits used in mutipolar expansion.

        Parameters
        ----------
        pcell : pheasy.atoms
            Unit cell structure.
        symops : dict
            Symmetry operations of space group.
            It has the keys "translations" and "rotations",
            and the values has the length of total number of symmetry.
        max_order : int
            highest order in generating multipole space.
        eps : float
            Numerical tolerance for symmetry operations.

        Returns
        -------
        MultipoleSpace
            An instance of class MultipoleSpace containing all the symmetry-
            distinct sites and their orbits for mulipolar expansion of charge
            density response.

        """
        nsym = pcell.symops["translations"].shape[0]
        atoms = pcell.scaled_positions
        types = pcell.get_atomic_numbers()
        symbols = pcell.get_chemical_symbols()
        if "equivalent_atoms" in pcell.symops:
            equivalent_atoms = pcell.symops["equivalent_atoms"]
            _, distict_atom_index = np.unique(equivalent_atoms, return_index=True)
        else:
            distict_atom_index = np.arange(
                pcell.get_global_number_of_atoms(), dtype=int
            )
        print("Atoms to be calculated:")

        MS = []
        orbit_space = set()
        for rep_idx in distict_atom_index:
            if rep_idx not in orbit_space:
                clus_orbit = deque()
                isotropy = deque()
                repeated = set()
                rep_clus = Site(
                    atoms[rep_idx],
                    rep_idx,
                    types[rep_idx],
                    symbols[rep_idx],
                    lrep=True,
                )
                for s in range(nsym):
                    trans = pcell.symops["translations"][s]
                    # if np.linalg.norm(trans) > eps: continue
                    rotmat = pcell.symops["rotations"][s]
                    hid_atom = np.dot(rotmat, atoms[rep_idx]) + trans
                    hid_idx = np.where(
                        np.linalg.norm(
                            atoms - hid_atom - np.rint(atoms - hid_atom),
                            axis=1,
                        )
                        < eps
                    )[0][0]
                    if rep_idx == hid_idx:
                        crotmat = get_rotation_cartesian(rotmat, pcell.cell)
                        iso_clus = {
                            "rotmat": rotmat,
                            "crotmat": crotmat,
                            "trans": trans,
                        }
                        isotropy.append(iso_clus)
                    if hid_idx not in repeated:
                        repeated.add(hid_idx)
                        orbit_space.add(hid_idx)
                        hid_clus = Site(
                            hid_atom,
                            hid_idx,
                            types[rep_idx],
                            symbols[rep_idx],
                            rotmat=rotmat,
                            trans=trans,
                            lrep=False,
                        )
                        hid_clus.crotmat = get_rotation_cartesian(
                            hid_clus.rotmat, pcell.cell
                        )
                        clus_orbit.append(hid_clus)
                rep_clus.set_isotropy_group(list(isotropy))
                clus_orbit.appendleft(rep_clus)
                MS.append(list(clus_orbit))
                print(f"- Index {rep_idx} | {symbols[rep_idx]}")

        return MultipoleSpace(pcell.cell, pcell.symops, MS, max_order)

    def set_multipole_tensors(self, mp_mat, natom):
        """Set the multipole tensors as an attribute.

        Parameters
        ----------
        mp_mat : numpy.ndarray
            A flattened array of full multipole tensors for
            symmetry-distinctive atoms.
        natom : int
            The number of atoms in (primitive) unit cell.

        """
        clusters = self.get_multipole_space()
        ncomp = self.get_number_of_components_all_order()

        multipole = OrderedDict()
        for order in range(2, self.max_order + 1):
            shape = [natom] + [3] * order
            size = np.power(3, order)
            mp_mat_full = np.zeros(shape)
            istart = np.sum(ncomp[: order - 2], dtype=int)
            iend = np.sum(ncomp[: order - 1], dtype=int)
            mp_mat_order = mp_mat[istart:iend]
            for n, orbit in enumerate(clusters):
                for cluster in orbit[1:]:
                    ind = cluster.atom_index
                    Gamma = cluster.get_crotation_tensor(order).toarray()
                    tmp = Gamma.dot(mp_mat_order[n * size : (n + 1) * size])
                    mp_mat_full[ind, ...] = tmp.reshape([3] * order)
            multipole[order] = mp_mat_full
        self._natom = natom
        self._multipole = multipole

    def set_dielectric_tensors(self, epsil_order, epsil_mat):
        """Set the dielectric tensors as an attribute.

        Parameters
        ----------
        epsil_order : int
            The maximum order of dielectric tensor used to expand
            macroscopic dielectric function.
        epsil_mat : numpy.ndarray
            Dielectric tensors up to the specified order stored
            as a 1D array.

        """
        istart = 0
        iend = 0
        epsilon = OrderedDict()
        for order in range(2, epsil_order + 1, 2):
            iend += np.power(3, order)
            epsil_mat_order = epsil_mat[istart:iend]
            epsilon[order] = epsil_mat_order.reshape([3] * order)
            istart += np.power(3, order)
        self._epsilon = epsilon
        self._epsil_order = epsil_order

    def get_multipole_tensors(self, order=None, lrep=False):
        """Return the multipole tensors at the given order.

        Parameters
        ----------
        order : int
            The order of multipole tensors to be returned. If
            not set or None, the full multipole tensors up the
            maximum order used in expansion will be returned.
        lrep : bool
            If True, only the results of representative site will
            be returned.

        Returns
        -------
        numpy.ndarray
            The multipole tensor in the shape of (natom, 3, 3, ...) or
            (size, 3, 3, ...) if `lrep` is True at the given order.
            If `order` is None, the result of full multipole tensors
            will flattened into a 1D array.

        """
        if lrep:
            rep_sites = self.get_representative_sites()
            if order is not None:
                return self._multipole[order][rep_sites]
            else:
                for order in range(2, self._max_order + 1):
                    mp_arr = [_.flatten() for _ in self._multipole[order][rep_sites]]
                return np.hstack(mp_arr)
        else:
            if order is not None:
                return self._multipole[order]
            else:
                mp_arr = [_.flatten() for _ in self._multipole.values()]
                return np.hstack(mp_arr)

    def get_dielectric_tensors(self, order=None):
        """Return the dielectic tensors at the given order.

        Parameters
        ----------
        order : int
            The order of dielectric tensors to be returned.
            It should be even order. If not set or None, the
            full dielectric tensors up the maximum order used
            in expansion will be returned.


        Returns
        -------
        numpy.ndarray
            The dielectric tensors in shape (3, 3, ...) at
            the given order. If order is None, the result of full
            dielectric tensors will flattened into a 1D array.

        """
        if order is not None:
            return self._epsilon[order]
        else:
            epsil_arr = [_.flatten() for _ in self._epsilon.values()]
            return np.hstack(epsil_arr)

    def write_multipole_tensors(self):
        """Write multipole tensors into files.

        This method will automatically write into files all of calculated
        multipole tensors up to the maximum order used in expansion.

        """
        for order in range(2, self._max_order + 1):
            mp_mat_full = self._multipole[order]

            if order == 2:
                """Write Born effective charge tensors."""
                filename = "born_charge.fmt"
                with open(filename, "w") as fh:
                    for iat in range(self._natom):
                        fh.write(f"{iat+1:>5d}\n")
                        fh.write(
                            "".join(map(lambda x: f"{x:25.15f}", mp_mat_full[iat, 0]))
                            + "\n"
                        )
                        fh.write(
                            "".join(map(lambda x: f"{x:25.15f}", mp_mat_full[iat, 1]))
                            + "\n"
                        )
                        fh.write(
                            "".join(map(lambda x: f"{x:25.15f}", mp_mat_full[iat, 2]))
                            + "\n"
                        )

            elif order == 3:
                """Write quadrupole tensors."""
                filename = "quadrupole.fmt"
                with open(filename, "w") as fh:
                    fh.write(f"{'atom':^8s}{'dir':^5s}")
                    fh.write(f"{'Qxx':^34s}{'Qyy':^16s}{'Qzz':^34s}")
                    fh.write(f"{'Qyz':^16s}{'Qxz':^34s}{'Qxy':^16s}\n")
                    for iat in range(self._natom):
                        for idir in range(3):
                            fh.write(f"{iat+1:^8d}{idir+1:^5d}")
                            fh.write(f"{mp_mat_full[iat,idir,0,0]:>25.15f}")
                            fh.write(f"{mp_mat_full[iat,idir,1,1]:>25.15f}")
                            fh.write(f"{mp_mat_full[iat,idir,2,2]:>25.15f}")
                            fh.write(f"{mp_mat_full[iat,idir,1,2]:>25.15f}")
                            fh.write(f"{mp_mat_full[iat,idir,0,2]:>25.15f}")
                            fh.write(f"{mp_mat_full[iat,idir,0,1]:>25.15f}\n")

            elif order == 4:
                """Write octupole tensors."""
                filename = "octupole.fmt"
                with open(filename, "w") as fh:
                    for iat in range(self._natom):
                        fh.write(f"{iat+1:>3d}\n")
                        mp_part = mp_mat_full[iat, ...].flatten()
                        for i, ind in enumerate(itertools.product([1, 2, 3], repeat=4)):
                            fh.write(f"{ind[0]:5d}{ind[1]:5d}{ind[2]:5d}{ind[3]:5d}")
                            fh.write(f"{mp_part[i]:25.15f}\n")

            elif order == 5:
                """Write hexadecapole tensors."""
                filename = "hexadecapole.fmt"
                with open(filename, "w") as fh:
                    for iat in range(self._natom):
                        fh.write(f"{iat+1:>3d}\n")
                        mp_part = mp_mat_full[iat, ...].flatten()
                        for i, ind in enumerate(itertools.product([1, 2, 3], repeat=5)):
                            fh.write(
                                f"{ind[0]:5d}{ind[1]:5d}{ind[2]:5d}{ind[3]:5d}{ind[4]:5d}"
                            )
                            fh.write(f"{mp_part[i]:25.15f}\n")

            elif order == 6:
                """Write triacontadipole tensors."""
                filename = "triacontadipole.fmt"
                with open(filename, "w") as fh:
                    for iat in range(self._natom):
                        fh.write(f"{iat+1:>3d}\n")
                        mp_part = mp_mat_full[iat, ...].flatten()
                        for i, ind in enumerate(itertools.product([1, 2, 3], repeat=6)):
                            fh.write(f"{ind[0]:5d}{ind[1]:5d}{ind[2]:5d}")
                            fh.write(f"{ind[3]:5d}{ind[4]:5d}{ind[5]:5d}")
                            fh.write(f"{mp_part[i]:25.15f}\n")

    def read_multipole_tensors(self, natom, max_order, non_polar=False):
        """Read multipole tensors into files.

        Parameters
        ----------
        natom : int
            The number of atoms in (primitive) unit cell.
        max_order : int
            Read multipole tensors up to this order.
        non_polar : bool
            If True, Born effective charges will be set to zero.

        Returns
        -------
        numpy.ndarray
            An 1D array of multipole tensors for all of symmetry-
            distinct atoms in (primitive) unit cell up to the
            specified max_order.

        """
        atom_index = self.get_representative_sites()
        ncomp = self.get_number_of_components_all_order()

        mp_arr = np.zeros(np.sum(ncomp[: max_order - 1], dtype=int))
        for order in range(2, max_order + 1):
            mp_arr_order = np.zeros(ncomp[order - 2])
            size = np.power(3, order)
            istart = np.sum(ncomp[: order - 2], dtype=int)
            iend = np.sum(ncomp[: order - 1], dtype=int)

            if order == 2:
                """Read Born effective charge tensors."""
                if not non_polar:
                    with open("born_charge.fmt") as fh:
                        for iat in range(natom):
                            fh.readline()
                            mp_part = np.zeros([3] * order)
                            for idir in range(3):
                                line = fh.readline().split()
                                mp_part[idir, :] = np.array([float(_) for _ in line])
                            if iat in atom_index:
                                ind = atom_index.index(iat)
                                mp_arr_order[ind * size : (ind + 1) * size] = (
                                    mp_part.flatten()
                                )
            elif order == 3:
                """Read quadrupole tensors."""
                try:
                    with open("quadrupole.fmt") as fh:
                        fh.readline()  # skip comment line
                        for iat in range(natom):
                            mp_part = np.zeros([3] * order)
                            for idir in range(3):
                                line = fh.readline().split()
                                mp_part[idir, 0, 0] = float(line[2])
                                mp_part[idir, 1, 1] = float(line[3])
                                mp_part[idir, 2, 2] = float(line[4])
                                mp_part[idir, 1, 2] = float(line[5])
                                mp_part[idir, 0, 2] = float(line[6])
                                mp_part[idir, 0, 1] = float(line[7])
                            if iat in atom_index:
                                ind = atom_index.index(iat)
                                mp_arr_order[ind * size : (ind + 1) * size] = (
                                    mp_part.flatten()
                                )
                except FileNotFoundError:
                    print(
                        "File quadrupole.fmt not found, assuming vanishing"
                        + " quadrupole tensor by symmetry."
                    )
            elif order >= 4:
                """Read octupole, hexadecapole and triacontadipole tensors."""
                if order == 4:
                    filename = "octupole.fmt"
                elif order == 5:
                    filename = "hexadecapole.fmt"
                elif order == 6:
                    filename = "triacontadipole.fmt"
                with open(filename) as fh:
                    for iat in range(natom):
                        fh.readline()
                        mp_part = np.zeros([3] * order)
                        for _, item in enumerate(
                            itertools.product([0, 1, 2], repeat=order)
                        ):
                            mp_part[item] = float(fh.readline().split()[-1])
                        if iat in atom_index:
                            ind = atom_index.index(iat)
                            mp_arr_order[ind * size : (ind + 1) * size] = (
                                mp_part.flatten()
                            )

            mp_arr[istart:iend] = mp_arr_order

        return mp_arr

    def write_dielectric_tensors(self, filename="epsilon.fmt"):
        """Write dielectric tensors into files.

        This method will automatically write into files all of calculated
        dielectric tensors up to the maximum order used in expansion.

        Parameters
        ----------
        filename : str
            The filename where dielectric tensors are written.

        """
        with open(filename, "w") as fh:
            for order in range(2, self._epsil_order + 1, 2):
                epsil_mat_order = self._epsilon[order]
                if order == 2:
                    fh.write(f"{order:>3d}\n")
                    fh.write(
                        "".join(map(lambda x: f"{x:25.15f}", epsil_mat_order[0])) + "\n"
                    )
                    fh.write(
                        "".join(map(lambda x: f"{x:25.15f}", epsil_mat_order[1])) + "\n"
                    )
                    fh.write(
                        "".join(map(lambda x: f"{x:25.15f}", epsil_mat_order[2])) + "\n"
                    )
                else:
                    epsil_mat_order = epsil_mat_order.flatten()
                    fh.write(f"{order:>3d}\n")
                    for i, ind in enumerate(itertools.product([1, 2, 3], repeat=4)):
                        fh.write(f"{ind[0]:5d}{ind[1]:5d}{ind[2]:5d}{ind[3]:5d}")
                        fh.write(f"{epsil_mat_order[i]:25.15f}\n")

    @staticmethod
    def read_dielectric_tensors(epsil_order, filename="epsilon.fmt"):
        """Read dielectric tensors from files.

        Parameters
        ----------
        epsil_order : int
            The maximum order of dielectric tensor used to expand
            macroscopic dielectric function.
        filename : str
            The filename where dielectric tensors are written.

        Returns
        -------
        numpy.ndarray
            Dielectric tensors up to the specified order stored
            as a 1D array.

        """
        with open(filename, "r") as fh:
            epsil_arr = []
            for order in range(2, epsil_order + 1, 2):
                ord = int(fh.readline().split()[0])
                assert order == ord
                epsil_part = np.zeros([3] * order)
                if order == 2:
                    for idir in range(3):
                        line = fh.readline().split()
                        epsil_part[idir, :] = np.array([float(_) for _ in line])
                else:
                    for _, item in enumerate(
                        itertools.product([0, 1, 2], repeat=order)
                    ):
                        epsil_part[item] = float(fh.readline().split()[-1])
                epsil_arr.append(epsil_part.flatten())
        epsil_arr = np.hstack(epsil_arr)

        return epsil_arr

    def __len__(self):
        """Return the number of representative sites in this space."""
        return len(self._MS)

    def __repr__(self):
        """Return multipole space information."""
        idx = []
        symbol = []
        for item in self._MS:
            idx.append(item[0].atom_index)
            symbol.append(item[0].atom_symbol)
        txt = "MultipoleSpace(max_order: {!r}, atom_index: {!r}, atom_symbol: {!r})"
        return txt.format(self._max_order, idx, symbol)


class MultipoleSymmetry(object):
    """Class for building symmetry constraints for dynamical multipoles.

    The method 'impose_multipole_symmetry' should be called to impose
    all kinds of symmetry constraints on multipole tensors, including
    space group (isotropy) symmetry, permutation symmetry, acoustic sum
    rules for charge neutrality condition (i.e. on Born effective charges).
    If the system has vanishing Born effective charge tensor, then the
    charge neutrality condition must be enforced on dynamical quadrupole
    tensors.

    """

    """Global configuration for MultipoleSymmetry."""
    NullSpaceMultipoleFile = "ns_mp.npz"
    NullSpaceEpsilonFile = "ns_epsil.npz"

    def __init__(
        self,
        pcell,
        cry_basis=False,
        epsil_order=None,
        non_polar=False,
        eps=1e-4,
    ):
        """Initialize with essential information for symmtery constraints.

        Parameters
        ----------
        pcell : pheasy.Atoms
            Primitive unit cell.
        cry_basis : bool
            True to deal with rotation matrix in crystal coordinate,
            False for Cartesian coordinate.
        epsil_order : int
            The order of dielectric tensor used to expand macroscopic
            dielectric function. If None, do not construct the symmetry
            for dielectric tensor.
        non_polar : bool
            In non-polar solids, Born effective charges vanish. If True,
            remove the null space for Born effectives.
        eps : float
            Numeric tolerance used in null space construction.

        """
        self._eps = eps
        self._pcell = pcell
        self._cry_basis = cry_basis
        self._epsil_order = epsil_order
        self._non_polar = non_polar

    @classmethod
    def run(cls, settings, pcell, multipole_space):
        """Run symmetry constraints for multipole tensors.

        This is a class method that creates a MultipoleSymmetry
        instance and should be only called by the main program.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.
        pcell : Pheasy.Atoms
            Primitive unit cell.
        multipole_space : MultipoleSpace
            An instance that contains all symmetry-distinct sites
            followed by their orbits.

        Returns
        -------
        Tuple
            A tuple of null space for multipole and dielectric tensors.

        """
        NS_full = None
        NS_epsil = None

        if settings.MP_EXPAND:
            print("Starting to symmetrize multipole tensors.")
            if settings.CRY_BASIS:
                print("Symmmetry constraints are imposed in crystal coordinate.")
            else:
                print("Symmmetry constraints are imposed in Cartesian coordinate.")
            multipole_symmetry = cls(
                pcell,
                settings.CRY_BASIS,
                epsil_order=settings.EPSIL_ORDER,
                non_polar=settings.NON_POLAR,
                eps=settings.EPS,
            )
            NS_full = multipole_symmetry.impose_multipole_symmtery(multipole_space)
            sparse.save_npz(cls.NullSpaceMultipoleFile, NS_full)

            if settings.EPSIL_ORDER is not None:
                print("Symmmetrization of dielectric tensors is performed.")
                NS_epsil = multipole_symmetry.impose_dielectric_symmetry(
                    multipole_space
                )
                sparse.save_npz(cls.NullSpaceEpsilonFile, NS_epsil)
        else:
            print("Reconstructing null space of multipole symmetry from file.")
            NS_full = cls.construct_null_space_restart(
                settings.MP_ORDER, settings.NON_POLAR, settings.FIX_ORDER
            )
            if settings.EPSIL_ORDER is not None:
                print("Reconstructing null space of dielectric symmetry from file.")
                NS_epsil = cls.construct_null_space_dielectric_restart(
                    settings.EPSIL_ORDER, settings.FIX_EPSIL_ORDER
                )

        return NS_full, NS_epsil

    def get_number_of_free_parameters_tensor(self):
        """Return number of free components for multipoles of each order.

        Returns
        -------
        Dict
            It contains the keys of order with the values as a list
            of free parameters for multipole tensors of each order.

        """
        return self._ncomp_free_tensor

    def get_number_of_free_parameters_order(self):
        """Return number of free tensor components of each order.

        Returns
        -------
        Dict
            It contains the keys of order with the corresponding
            number of free components of multipole tensors of each
            order.

        """
        return self._ncomp_free_order

    def get_total_number_of_free_parameters(self):
        """Return total number of free components of multipole tensors."""
        return self._ncomp_free_tot

    @timeit("Construction of multipole symmetry")
    def impose_multipole_symmtery(self, multipole_space):
        """Build constraint matrix from all kinds of symmetry on mutipoles.

        Parameters
        ----------
        multipole_space : MultipoleSpace
            An instance that contains all symmetry-distinct representative
            sites for multipoles followed by their orbits.

        Returns
        -------
        scipy.sparse.coo_matrix
            Full null space in COO sparse matrix form.

        """
        max_order = multipole_space.max_order
        clusters = multipole_space.get_multipole_space()
        cons_mat_dict = {}
        null_space = {}
        null_space_list = []
        ncomp_free_tensor = {}
        ncomp_free_order = {}

        for order in range(2, max_order + 1):
            start_time = datetime.datetime.now()
            if order == 2:
                print("Symmetry constraints on Born effective charge tensors.")
            elif order == 3:
                print("Symmetry constraints on quadrupole tensors.")
            elif order == 4:
                print("Symmetry constraints on octupole tensors.")
            elif order == 5:
                print("Symmetry constraints on hexadecapole tensors.")
            elif order == 6:
                print("Symmetry constraints on triacontadipoles tensors.")
            else:
                print("Symmetry constraints on {order}-order multipole tensors.")
            cons_mat_dict[order] = deque()
            null_space[order] = deque()
            ncomp_free_tensor[order] = deque()
            ncomp_free_order[order] = deque()

            """Isotropy and permutation symmetry constraints."""
            print("- Imposing space group symmetry.")
            if order > 2:
                print("- Imposing permutation symmetry.")
            for orbit in clusters:
                cons_list = orbit[0].build_isotropy_symmetry_constraints(
                    order, self._cry_basis
                )
                if order > 2:
                    cons_list += orbit[0].build_permutation_symmetry_constraints(order)
                ns_mat = np.eye(3**order)
                for cons_mat in cons_list:
                    ns_mat_tmp = null_space_dense(cons_mat.dot(ns_mat), self._eps)
                    if ns_mat_tmp.shape == (0,):
                        ns_mat = ns_mat_tmp
                        break
                    else:
                        ns_mat = ns_mat.dot(ns_mat_tmp)
                        ns_mat = normalize(ns_mat, axis=0)
                ns_mat = sparse.coo_matrix(ns_mat)
                null_space[order].append(ns_mat)
                ncomp_free_tensor[order].append(ns_mat.shape[1])
                cons_mat_dict[order].append(cons_mat)

            """Charge neutrality condition."""
            if order == 2:
                asr = True
            elif order == 3 and ncomp_free_order[2] == 0:
                asr = True
            else:
                asr = False
            if asr:
                print("- Imposing charge neutrality condition.")
                cons_asr = self.build_charge_neutrality_condition(clusters, order)
                cons_mat_dict[order].append(cons_asr)

            """Construct whole null space matrix of each order."""
            print("- Calculating null space.")
            ns_mat = block_diag_sparse(null_space[order], order)
            if asr:
                ns_mat = ns_mat.toarray()
                ns_mat_tmp = null_space_dense(cons_asr.dot(ns_mat), self._eps)
                if ns_mat_tmp.shape == (0,):
                    ns_mat = ns_mat_tmp
                else:
                    ns_mat = ns_mat.dot(ns_mat_tmp)
                    ns_mat = normalize(ns_mat, axis=0)
                ns_mat = sparse.coo_matrix(ns_mat)
            ncomp_free_order[order] = ns_mat.shape[1]
            null_space_list.append(ns_mat)

            end_time = datetime.datetime.now()
            total_time = end_time - start_time
            print("- Time elapsed: {}.".format(total_time))

        # TODO: cry_basis
        if self._non_polar:
            ns_mat_full = sparse.block_diag(null_space_list[1:])
        else:
            ns_mat_full = sparse.block_diag(null_space_list)
        self._constraints = cons_mat_dict
        self._ncomp_free_tensor = ncomp_free_tensor
        self._ncomp_free_order = ncomp_free_order
        self._ncomp_free_tot = ns_mat_full.shape[1]

        print("Summary of multipole tensors:")
        for order in range(2, max_order + 1):
            size = np.power(3, order)
            if order == 2:
                print("- Born effective charges")
            elif order == 3:
                print("- Quadrupoles")
            elif order == 4:
                print("- Octupoles")
            elif order == 5:
                print("- Hexadecapoles")
            elif order == 6:
                print("- triacontadipoles")
            else:
                print(f"- {order}-order multipoles")
            print(
                f"  Size: {size * len(clusters):<3d} |"
                + f" Free: {ncomp_free_order[order]}"
            )
            for i, orbit in enumerate(clusters):
                idx = orbit[0].atom_index
                symbol = orbit[0].atom_symbol
                print(
                    f"  Index {idx:<3d} | {symbol:<3s} | Free: {ncomp_free_tensor[order][i]}"
                )
            sparse.save_npz(f"ns_mp{order}", null_space_list[order - 2])
        print("- Total number of free components: {}".format(self._ncomp_free_tot))

        return ns_mat_full

    @timeit("Construction of dielectric symmetry")
    def impose_dielectric_symmetry(self, multipole_space):
        """Build constraint matrix and null space for dielectric tensors.

        Parameters
        ----------
        multipole_space : MultipoleSpace
            An instance that contains all symmetry-distinct representative
            sites for multipoles followed by their orbits.

        Returns
        -------
        scipy.sparse.coo_matrix
            Full null space in COO sparse matrix form for dielectric tensors.

        """
        print("Symmetrizing dielectric tensors up to " + f"{self._epsil_order}-order.")
        symops = multipole_space.get_symmetry()
        cell = multipole_space.cell
        epsil_ns_list = self.symmetrize_dielectric_tensors(symops, cell)
        ns_mat_full = sparse.block_diag(epsil_ns_list)

        ncomp_free_order_epsil = {}
        for n, order in enumerate(range(2, self._epsil_order + 1, 2)):
            size = np.power(3, order)
            ncomp_free_order_epsil[order] = epsil_ns_list[n].shape[1]
            if order == 2:
                print("- Dielectric permittivity tensor")
            if order == 4:
                print("- Dielectric dispersion tensor")
            print(
                "  Size: {:<3d} | Free: {:<3d}".format(
                    size, ncomp_free_order_epsil[order]
                )
            )
            sparse.save_npz(f"ns_epsil{order}", epsil_ns_list[n])
        self._ncomp_free_order_epsil = ncomp_free_order_epsil
        self._ncomp_free_tot_epsil = np.sum(list(ncomp_free_order_epsil.values()))
        print(
            "- Total number of free components: {}".format(self._ncomp_free_tot_epsil)
        )

        return ns_mat_full

    def build_charge_neutrality_condition(self, clusters, order):
        """Construct symmetry constraints by charge neutrality.

        It only applies to Born effective charge tensors or
        dynamical quadrupole tensor. If the system has vanishing
        Born effective charge tensor, then the charge neutrality
        condition must be enforced on quadrupole tensors.

        Parameters
        ----------
        clusters : dict
            Cluster space for multipole expansion.
        order : int
            The order of clusters that symmetry constraints apply to.

        Returns
        -------
        scipy.sparse.coo_matrix
            Constraint matrix of acoustic sum rule for charge
            neutrality in COO sparse matrix form.

        """
        if self._cry_basis:
            rot_fmt = "rotmat"
        else:
            rot_fmt = "crotmat"
        size = np.power(3, order)

        cons_mat_tmp = deque()
        for orbit in clusters:
            block = sparse.coo_matrix((size, size))
            for cluster in orbit[1:]:
                Gamma = kron_product(getattr(cluster, rot_fmt), order)
                block += Gamma
            cons_mat_tmp.append(block)
        cons_mat = sparse.hstack(cons_mat_tmp)

        return cons_mat

    def symmetrize_dielectric_tensors(self, symops, cell):
        """Find independent components of dielectric tensors.

        Parameters
        ----------
        symops : Dict
            Symmetry operations of space group. It has the keys
            "translations" and "rotations", and the values has
            the length of total number of symmetry.
        cell : numpy.ndarray
            A (3, 3) matrix for lattice vectors of (primitive)
            unit cell.

        Returns
        -------
        List
            A list of null space for dielectric tensors up to the order
            specified by epsil_order.

        """
        print("- Imposing point group and permutation symmetry.")
        nsym = symops["translations"].shape[0]
        cons_mat_dict = dict()
        ns_mat_list = []
        for order in range(2, self._epsil_order + 1, 2):
            size = np.power(3, order)
            cons_mat_list = deque()

            """Space group or point group symmetry"""
            for s in range(nsym):
                # TODO: skip those operations with fractional translation?
                trans = symops["translations"][s]
                if np.linalg.norm(trans) > self._eps:
                    continue
                rotmat = symops["rotations"][s]
                if not self._cry_basis:
                    rotmat = get_rotation_cartesian(rotmat, cell)
                Gamma = kron_product(rotmat, order)
                cons_mat_list.append(Gamma - sparse.identity(size))

            """Permutation symmetry"""
            permuted_list = []
            for i in range(order - 1):
                for j in range(i + 1, order):
                    ind = list(range(order))
                    ind[i], ind[j] = ind[j], ind[i]
                    permuted_list.append(ind)
            ind = list(range(order))
            for item in permuted_list:
                Rmat = get_permutation_tensor(ind, item).reshape((size, size))
                cons_mat = sparse.identity(size) - sparse.coo_matrix(Rmat)
                cons_mat_list.append(cons_mat)

            ns_mat = np.eye(size)
            cons_mat_dict[order] = cons_mat_list
            for cons_mat in cons_mat_list:
                ns_mat_tmp = null_space_dense(cons_mat.dot(ns_mat), self._eps)
                if ns_mat_tmp.shape == (0,):
                    ns_mat = ns_mat_tmp
                    break
                else:
                    ns_mat = ns_mat.dot(ns_mat_tmp)
                    ns_mat = normalize(ns_mat, axis=0)
            ns_mat = sparse.coo_matrix(ns_mat)
            ns_mat_list.append(ns_mat)

        return ns_mat_list

    @staticmethod
    def construct_null_space_restart(max_order, non_polar=False, fix_order=None):
        """Construct full null space for multipoles from a restart.

        Through this method, it will build full null space matrix
        from the saved npz files of multipoles at each order.

        Parameters
        ----------
        max_order : int
            Highest order considered in multipole expansion.
        non_polar : bool
            If True, null space for Born effective charges is
            removed.
        fix_order : int
            Mulipole tensors up to this order are fixed during
            fitting procedure and thus those null space blocks
            are removed.

        Returns
        -------
        scipy.sparse.coo_matrix
            Full null space in COO sparse matrix form. It consists
            of null space for all multipoles.

        """
        if fix_order is not None:
            min_order = fix_order
        else:
            min_order = 1

        null_space_list = []
        print("Summary of multipole tensors:")
        for order in range(2, max_order + 1):
            null_space_order = sparse.load_npz(f"ns_mp{order}.npz")
            if order == 2:
                print("- Born effective charges")
            elif order == 3:
                print("- Quadrupoles")
            elif order == 4:
                print("- Octupoles")
            elif order == 5:
                print("- Hexadecapoles")
            elif order == 6:
                print("- triacontadipoles")
            else:
                print(f"- {order}-order multipoles")
            print(
                f"  Size: {null_space_order.shape[0]:<3d} |"
                + f" Free: {null_space_order.shape[1]:<3d}"
            )
            if non_polar and order == 2:
                continue
            if order > min_order:
                null_space_list.append(null_space_order)
        ns_mat_full = sparse.block_diag(null_space_list)
        print(f"- Total number of free components: {ns_mat_full.shape[1]}")

        return ns_mat_full

    @staticmethod
    def construct_null_space_dielectric_restart(epsil_order, fix_order=None):
        """Construct full null space for dielectric tensors from a restart.

        Through this method, it will build full null space matrix
        from the saved npz files of multipoles at each order.

        Parameters
        ----------
        epsil_order : int
            Highest order considered in expanding dielectric function.
        fix_order : int
            Dielectric tensors up to this order are fixed during
            fitting procedure and thus those null space blocks
            are removed.

        Returns
        -------
        scipy.sparse.coo_matrix
            Full null space in COO sparse matrix form. It consists
            of null space for dielectric tensors at each order.

        """
        if fix_order is not None:
            min_order = fix_order
        else:
            min_order = 1

        null_space_list = []
        print("Summary of dielectric tensors:")
        for order in range(2, epsil_order + 1, 2):
            null_space_order = sparse.load_npz(f"ns_epsil{order}.npz")
            if order == 2:
                print("- Dielectric permittivity tensor")
            elif order == 4:
                print("- Dielectric dispersion tensor")
            print(
                f"  Size: {null_space_order.shape[0]:<3d} |"
                + f" Free: {null_space_order.shape[1]:<3d}"
            )
            if order > min_order:
                null_space_list.append(null_space_order)
        ns_mat_full = sparse.block_diag(null_space_list)
        print(f"- Total number of free components: {ns_mat_full.shape[1]}")

        return ns_mat_full

    def write(self, filename="multipole_symmetry.pkl"):
        """Write MultipoleSymmetry instance into pickle file.

        Parameters
        ----------
        filename : str
            Filename to save multipole symmetry constraints.
            By default "multipole_symmetry.pkl".

        """
        with open(filename, "wb") as fd:
            pickle.dump(self, fd)

    @staticmethod
    def read(filename="multipole_symmetry.pkl"):
        """Read and create MultipoleSymmetry instance from pickle file.

        Parameters
        ----------
        filename : str
            Filename for MultipoleSymmetry instance to read from.

        Returns
        -------
        MultipoleSymmetry instance

        """
        with open(filename, "rb") as fd:
            return pickle.load(fd)


class MultipoleConstructor(object):
    """Main driver for constructing multipole tensors."""

    """Global configurations for multipole expansion."""
    QPointsFile = "qpoints.dat"
    SensingMatrixQ1File = "smq1_prime.npz"
    SensingMatrixQ2File = "smq2_prime.npz"
    ChargeDensityRespnseFile = "drho.npz"
    EffectiveChargeFile = "zeff.npz"
    EpsilonFunctionFile = "epsil.npz"
    DielectricTestFile = "dielectric_test.npz"
    DielectricPredictFile = "dielectric_pred.npz"
    PhFilePattern = "ph.in.{:02d}"
    DrhoFilePattern = "drho.dat.{:02d}"
    ZeffFilePattern = "effcharge1.dat.{:02d}"
    DrhodvFilePattern = "drhodv.dat.{:02d}"

    def __init__(self, settings):
        """Initialization function.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.

        """
        self._read_primitive_cell(settings)

        """Run multipole expansion."""
        self.multipole_space = MultipoleSpace.run(settings, self.pcell)

        """Enforce multipole and dielectric symmetry."""
        self.null_space, self.null_space_epsil = MultipoleSymmetry.run(
            settings, self.pcell, self.multipole_space
        )

        if settings.EPSIL_ORDER is not None:
            self._epsil_order = settings.EPSIL_ORDER
        else:
            self._epsil_order = None

    def _read_primitive_cell(self, settings):
        """Read primitive unit cell and analyze its symmetry.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.

        """
        self.pcell = Atoms.read(settings.CELL_FILENAME, alat=settings.ALAT)
        self.natom = self.pcell.get_global_number_of_atoms()

        """Print crystal system information."""
        space_group = self.pcell.get_space_group(symprec=settings.SYMPREC)
        symops = self.pcell.get_symmetry(symprec=settings.SYMPREC)
        nsym = self.pcell.get_number_of_symmetries()
        if settings.WRITE_SYMOPS:
            write_symops(symops, "cell.symops")
        print("System: %s" % self.pcell.get_chemical_formula())
        print(f"Space group: {space_group}, {nsym} symmetry operations found.")

    def run_sensing_matrix(self, settings):
        """A wrapper function for constructing sensing matrix of q perturbations.

        This method should be called in the main function of multipole exapnsion
        to either build sensing matrix from scratch or read it from file.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.

        """
        if settings.SENSING_MAT:
            """Create sensing matrix from scratch."""
            print("Starting to construct sensing matrix of multipole expansion.")

            if settings.Q_FILE:
                print("Reading q perturbations from file.")
                q_vecs = np.loadtxt(self.QPointsFile)
            else:
                rcell = self.pcell.reciprocal()
                if settings.MESH is not None:
                    print("Creating q perturbations on a mesh grid.")
                    alat_ang = self.pcell.get_lattice_constant(unit="Angstrom")
                    q_vecs = wavevector_perturb_uniform(
                        settings.MESH, settings.MESH_STEP
                    )
                    if settings.IR_Q:
                        print("Finding irreducible q points.")
                        symops = self.pcell.get_symmetry()
                        q_vecs = get_irreducible_qpoints(q_vecs, symops)
                    q_vecs = q_vecs.dot(rcell) * alat_ang
                else:
                    if settings.Q_DIR is not None:
                        print("Creating q perturbations along a specific direction.")
                        q_dir = rcell.T.dot(settings.Q_DIR)
                        q_dir = np.where(abs(q_dir) < settings.EPS, 0.0, q_dir)
                    else:
                        q_dir = None
                    q_vecs = wavevector_perturb_simple(
                        settings.NQPTS,
                        settings.Q_MAX,
                        settings.Q_MIN,
                        q_dir,
                        settings.RAND_SEED,
                    )
                np.savetxt(self.QPointsFile, q_vecs)

            self._compute_sensing_matrix(q_vecs)
            np.savez_compressed(self.SensingMatrixQ1File, mat=self._SMQ1_prime)
            if self._epsil_order is not None:
                np.savez_compressed(self.SensingMatrixQ2File, mat=self._SMQ2_prime)
        else:
            """Read sensing matrix from file."""
            if settings.FIT or settings.PREDICT:
                print("Reconstructing sensing matrix of multipole expansion from file.")
                self._SMQ1_prime = np.load(self.SensingMatrixQ1File)["mat"]
                if self._epsil_order is not None:
                    self._SMQ2_prime = np.load(self.SensingMatrixQ2File)["mat"]
                if settings.NQPTS is not None:
                    self._nqpts = settings.NQPTS
                    self._SMQ1_prime = self._SMQ1_prime[
                        : 3 * self.natom * settings.NQPTS
                    ]
                    if self._epsil_order is not None:
                        self._SMQ2_prime = self._SMQ2_prime[: settings.NQPTS]
                else:
                    self._nqpts = self._SMQ1_prime.shape[0] // (3 * self.natom)

    @timeit("Construction of sensing matrix for multipole expansion")
    def _compute_sensing_matrix(self, q_vecs):
        """Calculate the sensing matrix for q perturbations.

        Parameters
        ----------
        q_vecs : numpy.ndarray
            A list of wave vectors for q perturbations.

        """
        nqpts = q_vecs.shape[0]
        alat = self.pcell.get_lattice_constant(unit="Bohr")

        sensing_mat_q1 = deque()
        sensing_mat_q2 = deque()
        for iq in range(nqpts):
            q_vec = q_vecs[iq]
            q_norm = np.linalg.norm(q_vec)
            filename = self.PhFilePattern.format(iq + 1)
            write_ph_input(q_vec, filename)

            print(f"- Generating q-point {iq + 1} of {nqpts}")
            print(
                f"{q_vec[0]:>10.6f} {q_vec[1]:>10.6f} {q_vec[2]:>10.6f}"
                + f" |q|={q_norm:>.6f}."
            )
            q_vec = q_vec * 2.0 * np.pi / alat
            sensing_mat = build_multipole_wavevector_matrix(
                self.multipole_space, self.natom, q_vec
            )
            sensing_mat_q1.append(sensing_mat)
            if self._epsil_order is not None:
                sensing_mat = build_dielectric_wavevector_matrix(
                    self._epsil_order, q_vec
                )
                sensing_mat_q2.append(sensing_mat)

        self._SMQ1_prime = np.vstack(sensing_mat_q1)
        if self._epsil_order is not None:
            self._SMQ2_prime = np.vstack(sensing_mat_q2)

    def fit_multipole_expansion(self, settings):
        """Fit multipole and dielectric expansions using a linear model.

        This method should be called in the main function of multipole
        expansion to fit the multipole and dielectric expansions. Currently,
        only the least-square fitting is supported.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.

        """
        self._epsilon = None
        self._mp_arr_fix = None
        self._epsil_arr_fix = None

        exits_drhodv = os.path.isfile(self.DrhodvFilePattern.format(1))
        if settings.EPSIL_ORDER is not None:
            if not exits_drhodv:
                raise FileNotFoundError("Dielectric response function not found.")
            print("Starting to fit multipole and dielectric tensors by a linear model.")
        else:
            print("Starting to fit multipole tensors by a linear model.")
        if exits_drhodv:
            self._read_dielectric_response(
                self._nqpts,
                is_screened=settings.SCREENED,
                epsil_kernel=settings.EPSIL_KERNEL,
            )
        self._read_charge_density_response(self._nqpts)
        self._preprocess_sensing_matrix(
            non_polar=settings.NON_POLAR,
            fix_order=settings.FIX_ORDER,
            fix_epsil_order=settings.FIX_EPSIL_ORDER,
        )

        """Perform linear regression."""
        self._fit_multipole_tensors()
        if self._epsil_order is not None:
            self._fit_dielectric_tensors()

    @timeit("Multipole tensor fitting")
    def _fit_multipole_tensors(self):
        """Perform linear regression to fit multipole tensors."""
        print("Fitting multipole tensors via the ordinary least-square.")

        optimizer = Optimizer()
        optimizer.fit(self._SMQ1, self._drho)
        fit_results = optimizer.results
        fit_metrics = optimizer.metrics

        print("Summary of multipole tensors fitting:")
        print(f"- Free multipole components: {fit_results['n_featrues']}")
        print(f"- Rank of coefficient matrix: {fit_metrics['rank']}")
        print(f"- RMSE: {fit_metrics['rmse']} |e| / bohr^3")
        print(f"- Mean relative error: {fit_metrics['mape']}")
        print(f"- Relative norm error: {fit_metrics['rne']}")
        print(f"- r2 score: {fit_metrics['r2_score']}")

        """Write multipole tensors into file."""
        print("Writing multipole tensors into files.")
        mp_arr = self.null_space.dot(fit_results["parameters"])
        if self._mp_arr_fix is not None:
            mp_arr = np.hstack([self._mp_arr_fix, mp_arr])
        self.multipole_space.set_multipole_tensors(mp_arr, self.natom)
        self.multipole_space.write_multipole_tensors()

    @timeit("Dielectric tensor fitting")
    def _fit_dielectric_tensors(self):
        """Perform linear regression to fit dielectric tensors."""
        print("Fitting dielectric tensors via the ordinary least-square.")

        optimizer = Optimizer()
        optimizer.fit(self._SMQ2, self._xi)
        fit_results = optimizer.results
        fit_metrics = optimizer.metrics

        print("Summary of dielectric tensors fitting:")
        print(f"- Free dielectric components: {fit_results['n_featrues']}")
        print(f"- Rank of coefficient matrix: {fit_metrics['rank']}")
        print(f"- RMSE: {fit_metrics['rmse']}")
        print(f"- Mean relative error: {fit_metrics['mape']}")
        print(f"- Relative norm error: {fit_metrics['rne']}")
        print(f"- r2 score: {fit_metrics['r2_score']}")

        """Write dielectric tensors into file."""
        print("Writing dielectric tensors into files.")
        epsil_arr = self.null_space_epsil.dot(fit_results["parameters"])
        if self._epsil_arr_fix is not None:
            epsil_arr = np.hstack([self._epsil_arr_fix, epsil_arr])
        self.multipole_space.set_dielectric_tensors(self._epsil_order, epsil_arr)
        self.multipole_space.write_dielectric_tensors()

    def _read_dielectric_response(
        self, nqpts, is_screened=False, epsil_kernel=1, save_to_file=True
    ):
        r"""Read dielectric response function from file.

        Parameters
        ----------
        nqpts : int
            Number of q points.
        is_screened : bool, optional
            If True, the macroscopic potential is included in QE-DFPT,
            i.e. drho being the total charge density response. This should
            set to True when 'lmacro=.FALSE.' set in DFPT, otherwise False.
        epsil_kernel : int, optional
            The kernel for dielectric screening function. The allowed values
            are 1 for $1+v(q)*\chi(q)$ and 2 for $1/epsilon(q)^{-1}$. This
            option is only available when the macroscopic potential is not
            removed, i.e. 'lmacro=.FALSE.' in DFPT, otherwise the kernel
            $1/epsilon(q)^{-1}$ is used. By default 1.
        save_to_file : bool, optional
            If True, save the dielectric response function to file. By default True.

        """
        print(f"Reading inverse dielectric function, {nqpts} q-points.")
        alat = self.pcell.get_lattice_constant(unit="Bohr")
        q_vecs = np.loadtxt(self.QPointsFile) * 2.0 * np.pi / alat
        q_norms = np.linalg.norm(q_vecs, axis=1)

        xi = np.zeros(nqpts)
        if is_screened:
            epsilon = np.zeros(nqpts, dtype=complex)
        else:
            epsilon = None
        for iq in range(nqpts):
            filename = self.DrhodvFilePattern.format(iq + 1)
            drhodv = read_inverse_dielectric_function(filename)
            if not is_screened:
                xi[iq] = drhodv[-1].real * q_norms[iq] ** 2
            else:
                if epsil_kernel == 1:
                    epsilon[iq] = drhodv[-1]
                else:
                    epsilon[iq] = drhodv[1]
                xi[iq] = 1.0 / epsilon[iq].real * q_norms[iq] ** 2
            print(f"- {filename} read successfully.")
        if save_to_file:
            np.savez_compressed(self.EpsilonFunctionFile, epsil=xi)

        self._xi = xi
        self._epsilon = epsilon

    def _read_charge_density_response(self, nqpts, save_to_file=True):
        """Read charge density response function from file.

        Parameters
        ----------
        nqpts : int
            Number of q points.
        save_to_file : bool, optional
            If True, save the charge density response function to file.
            By default True.

        """
        print(f"Reading charge density response, {nqpts} q-points.")
        volume = self.pcell.get_volume("Bohr")

        drho_cpl = np.zeros([nqpts, self.natom, 3], dtype=complex)
        for iq in range(nqpts):
            filename = self.DrhoFilePattern.format(iq + 1)
            drho_cpl[iq] = read_charge_density_response(self.natom, filename)
            if self._epsilon is not None:
                drho_cpl[iq] /= self._epsilon[iq]
            print(f"- {filename} read successfully.")
        drho_cpl = drho_cpl.flatten() * volume / 2**0.5
        if save_to_file:
            np.savez_compressed(self.ChargeDensityRespnseFile, drho=drho_cpl)

        self._drho = drho_cpl

    def _read_effective_charges(self, nqpts, save_to_file=True):
        """Read effective charges from file.

        Parameters
        ----------
        nqpts : int
            Number of q points.
        save_to_file : bool, optional
            If True, save the effective charges to file. By default True.

        """
        print(f"Reading effective charges, {nqpts} q-points.")

        zeff = np.zeros([nqpts, self.natom, 3], dtype=complex)
        for iq in range(nqpts):
            filename = self.ZeffFilePattern.format(iq + 1)
            zeff[iq] = read_charge_density_response(self.natom, filename)
            if self._epsilon is not None:
                zeff[iq] /= self._epsilon[iq]
            print(f"- {filename} read successfully.")
        if save_to_file:
            np.savez_compressed(self.EffectiveChargeFile, zeff=zeff)

        self._zeff = zeff

    def _preprocess_sensing_matrix(
        self, non_polar=False, fix_order=None, fix_epsil_order=None
    ):
        """Preprocess sensing matrix for fitting multipole tensors.

        This method should be called before performing a linear regression.
        It removes the parts of the fixed order of multipole and dielectric
        expansions from the sensing matrix.

        Parameters
        ----------
        non_polar : bool, optional
            Set to True if the system is non-polar where the Born effective
            charges are set to zero. By default False.
        fix_order : int, optional
            The order up to which the multipole expansion is fixed
            during the fitting. By default None.
        fix_epsil_order : int, optional
            The order up to which the dielectric expansion is fixed
            during the fitting. By default None.

        """
        if self._epsil_order is not None:
            if fix_epsil_order is not None:
                print(
                    "Preprocessing sensing matrix of multipole and dielectric expansions."
                )
                print(f"Dielectric tensors up to {fix_epsil_order}-order kept fixed.")
                nskip = 0
                for order in range(2, fix_epsil_order + 1, 2):
                    nskip += np.power(3, order)
                SMQ2_prime_fix = self._SMQ2_prime[:, :nskip]
                self._epsil_arr_fix = self.multipole_space.read_dielectric_tensors(
                    fix_epsil_order
                )
                self._xi -= SMQ2_prime_fix.dot(self._epsil_arr_fix)
                self._SMQ2_prime = self._SMQ2_prime[:, nskip:]
            self._SMQ2 = self._SMQ2_prime.dot(self.null_space_epsil.toarray())

        if fix_order is not None:
            print("Preprocessing sensing matrix of multipole expansion.")
            if non_polar:
                print(
                    "Non-polar system detected, "
                    + "assuming vanishing Born effective charges."
                )
            print(f"Multipole tensors up to {fix_order}-order kept fixed.")
            self._mp_arr_fix = self.multipole_space.read_multipole_tensors(
                self.natom, fix_order, non_polar
            )
            nskip = self.multipole_space.get_number_of_components_all_order()[
                : fix_order - 1
            ]
            nskip = np.sum(nskip, dtype=int)
            SMQ1_prime_fix = self._SMQ1_prime[:, :nskip]
            self._drho -= SMQ1_prime_fix.dot(self._mp_arr_fix)
            self._SMQ1_prime = self._SMQ1_prime[:, nskip:]
        elif non_polar:
            print(
                "Non-polar system detected, assuming vanishing Born effective charges."
            )
            nskip = self.multipole_space.get_number_of_components_all_order()[0]
            self._mp_arr_fix = np.zeros(nskip)
            self._SMQ1_prime = self._SMQ1_prime[:, nskip:]
        SMQ1_cpl = self._SMQ1_prime.dot(self.null_space.toarray())
        self._SMQ1 = np.vstack([SMQ1_cpl.real, SMQ1_cpl.imag])
        self._drho = np.hstack([self._drho.real, self._drho.imag])

    def predict_dielectic_properties(self, settings):
        """Predict multipole and dielectric responses.

        This method serves as a wrapper function for predicting charge
        density response, effective charges and dielectric function.

        Parameters
        ----------
        settings : namespace defined by argparse
            User settings.

        """
        q_vecs = np.loadtxt(self.QPointsFile)
        q_norms = np.linalg.norm(q_vecs, axis=1)

        if settings.EPSIL_ORDER is not None:
            """Load or read dielectric tensors."""
            if not hasattr(self.multipole_space, "_epsilon") or settings.MP_READ:
                print("Reading dielectric tensors from file.")
                epsil_arr = self.multipole_space.read_dielectric_tensors(
                    settings.EPSIL_ORDER
                )
                self.multipole_space.set_dielectric_tensors(
                    settings.EPSIL_ORDER, epsil_arr
                )

        """Load or read multipole tensors."""
        if not hasattr(self.multipole_space, "_multipole") or settings.MP_READ:
            print("Reading multipole tensors from file.")
            mp_arr = self.multipole_space.read_multipole_tensors(
                self.natom, settings.MP_ORDER
            )
            self.multipole_space.set_multipole_tensors(mp_arr, self.natom)

        """Read ab initio inverse macroscopic dielectric function,
           charge density response and effective charges."""
        exits_drhodv = os.path.isfile(self.DrhodvFilePattern.format(1))
        if exits_drhodv:
            self._read_dielectric_response(
                q_vecs.shape[0],
                is_screened=settings.SCREENED,
                epsil_kernel=settings.EPSIL_KERNEL,
                save_to_file=False,
            )
        self._read_charge_density_response(q_vecs.shape[0], save_to_file=False)
        self._read_effective_charges(q_vecs.shape[0], save_to_file=False)
        np.savez_compressed(
            self.DielectricTestFile,
            q=q_norms,
            drho=self._drho,
            zeff=self._zeff,
            xi=self._xi,
        )

        """Predict charge density response, effective charges and
           dielectric function from the fitted expansion."""
        pred_dict = self._predict_dielectric_response(q_vecs)
        np.savez_compressed(self.DielectricPredictFile, **pred_dict)

    def _predict_dielectric_response(self, qpts):
        """Predict dielectric response functions from the fitted expansion.

        TODO: implement the prediction using Einstein summation.

        Parameters
        ----------
        qpts : numpy.ndarray
            A list of wave vectors for q perturbations with the shape of (nqpts, 3).

        Returns
        -------
        dict
            A dictionary containing the predicted dielectric properties.

        """
        print("Predicting dielectric properties using the fitted expansion.")
        nqpts = qpts.shape[0]
        alat = self.pcell.get_lattice_constant(unit="Bohr")
        q_norms = np.linalg.norm(qpts, axis=1) * 2 * np.pi / alat
        volume = self.pcell.get_volume("Bohr")

        if self._epsil_order is not None:
            zeff_pred, xi_pred = self.multipole_space.predit(
                self.SMQ1_prime, self.SMQ2_prime
            )
            xi_pred /= q_norms**2
        else:
            zeff_pred = self.multipole_space.predit(self.SMQ1_prime)
            xi_pred = None
        zeff_pred = zeff_pred.reshape((nqpts, self.natom, 3))
        drho_pred = zeff_pred * np.sqrt(2) / volume
        zeff_pred /= q_norms[:, np.newaxis, np.newaxis] * -1j

        pred_dict = {
            "q": q_norms * alat / 2 / np.pi,
            "drho": drho_pred,
            "zeff": zeff_pred,
            "xi": xi_pred,
        }

        return pred_dict


def main():
    """Main function for running multipole expansion."""
    start_time = datetime.datetime.now()

    parser = InputParser()  # instantiate an input parser

    welcome(start_time)  # print welcome message

    """Parse user settings."""
    settings = parser.settings  # get settings from parser
    multipole = MultipoleConstructor(settings)  # instantiate a multipole constructor

    """Generate q perturbations and construct sensing matrix."""
    multipole.run_sensing_matrix(settings)

    """Fit multipole and dielectric expansions."""
    if settings.FIT:
        multipole.fit_multipole_expansion(settings)

    """Predict dielectric properties using fitted expansions."""
    if settings.PREDICT:
        multipole.predict_dielectic_properties(settings)

    """Finalize and estimate time cost"""
    goodbye(start_time)


if __name__ == "__main__":
    main()
