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
elec_positions = np.array([np.array([-50, 0, -200], dtype=float),
                           np.array([60, 126, 659], dtype=float),
                           np.array([-301, 39, 879], dtype=float)])


# elec_positions = [np.array([[50, 0, -200], ], dtype=float),
#                   np.array([[111, 131, 652], ], dtype=float),
#                   np.array([[-242, 7, 929], ], dtype=float),
#                   np.array([[-200, 0, 0], ], dtype=float),
#                   np.array([[50, 40, 120], ], dtype=float)]

current_amps = [1e4, -1e4]
# elec_positions = [np.array([[-242, 7, 929], ], dtype=float)]

# measure_coordinates = np.array([[0, 0, -200], [130, 131, 652],
#                                 [-242, 43, 929]])

measure_coordinates = np.array(
    [[-0, 0, - 200], [10, 126, 659], [-251, 39, 879]])
start = time.time()
# monophasic_pulse_params['stop_time'] = 20
extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params,
                      current_amps, measure_coordinates, 20)
end = time.time()
print(f'Time to execute {end - start} seconds')
