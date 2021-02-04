from main import ExternalPotentialSim
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
import time

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'Hallermann_ext_pot_test'

extPotSim = ExternalPotentialSim(cellsim_Hallermann_params)


# Test parameters
# current_amps = [2e4, -2e4, 1e4, -1e4, -8e3, 8e3, -7e3, 7e3, -5e3, 5e3]  # uA
elec_positions = np.array([[-50, 0, -200],
                           [60, 126, 659],
                           [-301, 39, 879]])


current_amps = [1e4, -1e4]

# measure_coordinates = np.array([[0, 0, -200], [130, 131, 652],
#                                 [-242, 43, 929]])

measure_coordinates = np.array(
    [[-0, 0, - 200], [10, 126, 659], [-251, 39, 879]])


def set_electrode_pos(measure_coordinates):
    """
    Sets electrodes at a given distance from the measurement coordinates
    Parameters:
    measure_coordinates
    Returns:
    """
    elec_positions = np.copy(measure_coordinates)
    for mc in range(elec_positions.shape[0]):

        elec_positions[mc][0] -= 50

    return elec_positions


ep = set_electrode_pos(measure_coordinates)
start = time.time()

for I in current_amps:

    for idx, pos in enumerate(elec_positions):

        extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params,
                              I, measure_coordinates, 20, pos, idx)
        # Plot Sim

end = time.time()
print(f'Time to execute {end - start} seconds')
