import numpy as np


def set_electrode_pos(measure_coordinates, shift1=-50, shift2=0):
    """
    Sets electrodes at a given distance from the measurement coordinates
    Parameters:
    measure_coordinates
    Returns:
        elec_positions, electrode positions near the compartments measurement
    """
    elec_positions = np.copy(measure_coordinates)

    for mc in range(elec_positions.shape[0]):

        elec_positions[mc][0] += shift1
        elec_positions[mc][2] += shift2

    return elec_positions
