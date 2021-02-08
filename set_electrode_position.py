import numpy as np


def set_electrode_pos(measure_coordinates):
    """
    Sets electrodes at a given distance from the measurement coordinates
    Parameters:
    measure_coordinates
    Returns:
    """
    print(measure_coordinates.shape)
    elec_positions = np.copy(measure_coordinates)
    print(elec_positions.shape)
    for mc in range(elec_positions.shape[0]):

        elec_positions[mc][0] -= 50

    return elec_positions
