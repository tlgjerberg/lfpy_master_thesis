from main import ExternalPotentialSim
from plotting import PlotSimulations
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
    print(measure_coordinates.shape)
    elec_positions = np.copy(measure_coordinates)
    print(elec_positions.shape)
    for mc in range(elec_positions.shape[0]):

        elec_positions[mc][0] -= 50

    return elec_positions


def run_hallermann(cell_models_folder, measure_coordinates, run_sim=False, plot_sim=False):

    start = time.time()
    extPotSim = ExternalPotentialSim(cellsim_Hallermann_params)

    if run_sim:
        elec_positions = set_electrode_pos(measure_coordinates)

        for I in current_amps:

            for idx, pos in enumerate(elec_positions):

                extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params,
                                      I, measure_coordinates, 20, pos, idx)
            # Plot Sim
    if plot_sim:
        cell_vmem = np.load(join(extPotSim.save_folder,
                                 'Hallermann_x_shift=0_z_rot=0_-10000.0mA_elec_pos=[ -50    0 -200]_vmem.npy'))
        cell_tvec = np.load(join(extPotSim.save_folder,
                                 'Hallermann_x_shift=0_z_rot=0_-10000.0mA_elec_pos=[ -50    0 -200]_tvec.npy'))
        plotSim = PlotSimulations(
            cellsim_Hallermann_params, cell_vmem, cell_tvec)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(cell, measure_coordinates)

    end = time.time()
    print(f'Time to execute {end - start} seconds')


run_hallermann(cell_models_folder, measure_coordinates, False, True)
