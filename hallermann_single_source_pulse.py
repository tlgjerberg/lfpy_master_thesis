from main import ExternalPotentialSim
from plotting import PlotSimulations
from set_electrode_position import set_electrode_pos
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
from mpi4py import MPI

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'Hallermann_ext_pot_test'


# Test parameters
# current_amps = [2e4, -2e4, 1e4, -1e4, -8e3, 8e3, -7e3, 7e3, -5e3, 5e3]  # uA
elec_positions = np.array([[-50, 0, -200],
                           [60, 126, 659],
                           [-301, 39, 879]])


current_amps = [1e4, -1e4]

measure_coordinates = np.array(
    [[-0, 0, - 200], [10, 126, 659], [-251, 39, 879]])


def run_hallermann(cell_models_folder, measure_coordinates, run_sim=False, plot_sim=False):

    start = time.time()
    extPotSim = ExternalPotentialSim(
        cellsim_Hallermann_params, monophasic_pulse_params)

    if run_sim:
        elec_positions = set_electrode_pos(measure_coordinates)

        extPotSim.run_ext_sim(cell_models_folder, I,
                              measure_coordinates, 20, pos)
        # Plot Sim
    if plot_sim:
        cell_vmem = np.load(join(extPotSim.save_folder,
                                 f'Hallermann_x_shift=0_z_rot=0_{I}mA_elec_pos={pos}_vmem.npy'))
        cell_tvec = np.load(join(extPotSim.save_folder,
                                 f'Hallermann_x_shift=0_z_rot=0_{I}mA_elec_pos={pos}_tvec.npy'))
        plotSim = PlotSimulations(
            cellsim_Hallermann_params, monophasic_pulse_params)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(cell, measure_coordinates, cell_vmem, cell_tvec)

    end = time.time()
    print(f'Time to execute {end - start} seconds')


task_idx = -1
for I in current_amps:
    for pos in elec_positions:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            # print('RANK: ', RANK, I, pos)
            continue
            # Lines below are only needed if NEURON needs to be reset every time
        pid = os.fork()
        if pid == 0:
            # Do work here!
            run_hallermann(cell_models_folder, measure_coordinates, True, True)
            print("RANK %d doing task %d" % (RANK, task_idx))
            os._exit(0)
        else:
            os.waitpid(pid, 0)
