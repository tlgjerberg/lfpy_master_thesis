from main import ExternalPotentialSim
from plotting import PlotSimulations
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_bisc_stick_params)
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
from mpi4py import MPI
import time
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


"""
Improve plotting

Create plot of Hallermann pyrmidal cell (4 figs?) showing electrode near soma,
dendrite and axon to demonstrate if activation is possible and at what current
amplitude.

Stimulate using current of 10 muA and go lower after (Histed et Al)

Create figure of axial current over distance along axon. Use simple stick model.

Improve plot_cellsim_alt for easy reading and page formatting.

Fix issue with splitting the workload properly between processes. Currently loop
divides work in a way that causes crashing

"""


def run_axon(cell_models_folder, measure_coordinates, I, pos, run_sim=False, plot_sim=False):

    start = time.time()
    extPotSim = ExternalPotentialSim(
        cellsim_bisc_stick_params, monophasic_pulse_params)

    # cell_name = extPotSim.return_sim_name()

    if run_sim:
        elec_positions = set_electrode_pos(measure_coordinates)

        extPotSim.run_ext_sim(cell_models_folder,
                              I, measure_coordinates, 20, pos)

    else:
        print('No simulation run!')
        # Plot Sim
    if plot_sim:

        #         cell_vmem = np.load(
        #             join(extPotSim.save_folder, cell_name + '.npy'))
        #         cell_tvec = np.load(
        #             join(extPotSim.save_folder, cell_name + '.npy'))

        cell_tvec = np.load(
            join(extPotSim.save_folder, f'axon_x_shift=0_z_rot=0_{I}mA_elec_pos={pos}_tvec' + '.npy'))
        cell_vmem = np.load(
            join(extPotSim.save_folder, f'axon_x_shift=0_z_rot=0_{I}mA_elec_pos={pos}_vmem' + '.npy'))
        plotSim = PlotSimulations(
            cellsim_bisc_stick_params, monophasic_pulse_params)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(
            cell, measure_coordinates, cell_vmem, cell_tvec)
    else:
        print('No plots generated!')

    end = time.time()
    print(f'Time to execute {end - start} seconds')


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
current_amps = [-1e4]  # uA
elec_positions = np.array([[0, 0, -50],
                           [0, 0, -100],
                           [0, 0, -200],
                           [0, 0, -400]], dtype=float)
cellsim_bisc_stick_params['save_folder_name'] = 'mpi_axon_test'
axon_measure_idxs = np.array(
    [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])
monophasic_pulse_params['stop_time'] = 200

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
            run_axon(cell_models_folder, axon_measure_idxs, I, pos, True, True)
            print("RANK %d doing task %d" % (RANK, task_idx))
            os._exit(0)
        else:
            os.waitpid(pid, 0)


# # Test parameters

# cell_models_folder = join(os.path.dirname(__file__), "cell_models")
# current_amps = [-1e4]  # uA
# positions = np.array([np.array([0, 0, -50], dtype=float),
#                       np.array([0, 0, -100], dtype=float),
#                       np.array([0, 0, -200], dtype=float),
#                       np.array([0, 0, -400], dtype=float),
#                       np.array([0, 0, -500], dtype=float)])

# cellsim_bisc_stick_params['save_folder_name'] = 'mpi_axon_test'
# axon_measure_idxs = np.array(
#     [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])
# monophasic_pulse_params['stop_time'] = 200

# extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)
#
# extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params, current_amps,
#                       positions, axon_measure_idxs, 200, passive=True)
# assert np.allclose(recvbuf, RANK)
