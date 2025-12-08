"""
Parallel computation helpers for multiprocessing support.

This module provides wrapper functions for parallelizing computation
across different atom types using Python's multiprocessing library.
"""

import numpy as np
from ComputeDynamicProperties import ComputeDynamicProperties


def compute_pdos_parallel(args):
    """
    Wrapper function for parallel PDOS computation.

    Args:
        args: Tuple containing (atom_type_index, atom_pos_dict, atom_vel_dict,
                               pos, vel, parameters, latt, omega)

    Returns:
        Tuple of (atom_type_index, vacf_non, vacf_output, pdos)
    """
    j, atom_pos_dict, atom_vel_dict, pos, vel, parameters, latt, omega = args

    if j == parameters['num_types'] + 1:
        pdos_cal = ComputeDynamicProperties(pos, vel, parameters, latt)
    else:
        pdos_cal = ComputeDynamicProperties(
            atom_pos_dict[j], atom_vel_dict[j], parameters, latt
        )

    vacf_non, vacf_output, pdos = pdos_cal.pdos(omega)
    return j, vacf_non, vacf_output, pdos


def compute_dynamic_structure_parallel(args):
    """
    Wrapper function for parallel dynamic structure computation.

    Args:
        args: Tuple containing (atom_type_index, atom_pos_dict, atom_vel_dict,
                               pos, vel, parameters, latt, omega, write_params)

    Returns:
        Tuple of (atom_type_index, fd_scale, Sv, S_intgr)
    """
    j, atom_pos_dict, atom_vel_dict, pos, vel, parameters, latt, omega = args

    if j == parameters['num_types'] + 1:
        dynamic_cal = ComputeDynamicProperties(pos, vel, parameters, latt)
    else:
        dynamic_cal = ComputeDynamicProperties(
            atom_pos_dict[j], atom_vel_dict[j], parameters, latt
        )

    fd_scale = dynamic_cal.calculate_intermediate_scattering()
    Sv = dynamic_cal.calculate_dynamic_structure(omega, fd_scale)
    S_intgr = dynamic_cal.Integrate_dynamic_structure(Sv)

    return j, fd_scale, Sv, S_intgr


def compute_van_hove_parallel(args):
    """
    Wrapper function for parallel Van Hove correlation computation.

    Args:
        args: Tuple containing (atom_type_index, atom_pos_dict, atom_vel_dict,
                               pos, vel, parameters, latt)

    Returns:
        Tuple of (atom_type_index, Gr_mean, shells, r)
    """
    j, atom_pos_dict, atom_vel_dict, pos, vel, parameters, latt = args

    if j == parameters['num_types'] + 1:
        van_hov_cal = ComputeDynamicProperties(pos, vel, parameters, latt)
    else:
        van_hov_cal = ComputeDynamicProperties(
            atom_pos_dict[j], atom_vel_dict[j], parameters, latt
        )

    Gr_mean, shells, r = van_hov_cal.calculate_van_hove_function()
    return j, Gr_mean, shells, r
