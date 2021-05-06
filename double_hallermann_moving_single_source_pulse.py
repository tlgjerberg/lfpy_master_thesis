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
cellsim_Hallermann_params['save_folder_name'] = 'Hallermann_double_test'

"""
Placeholder variables
"""
# elec_positions = np.array(
#     [[725, 0, -575], [800, 0, -575], [-800, 0, -575], [-725, 0, -575]])
elec_positions = np.array(
    [[725, 60, -548], [800, 60, -548], [-800, -50, -548], [-725, -50, -548]])

current_amps = [-1e4]

measure_coordinates = np.array([[680, 0, -548]])

# cell_rot = [np.pi, 0, np.pi / 2, np.pi / 3, np.pi / 4, 1.4]
# cell_rot = [np.pi / 2, (3. / 2) * np.pi]
cell_rot = [np.pi, 0]


def run_double_hallermann(cell_models_folder, measure_coords, I, pos, z, run_sim=False, plot_sim=False):
    print('RANK', RANK)
    print('z', z)
    test_pos = set_electrode_pos(measure_coords, 50)
    print(test_pos)
    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_Hallermann_params['z_rot'] = z

    extPotSim = ExternalPotentialSim(
        cellsim_Hallermann_params, monophasic_pulse_params)
    sim_name = extPotSim.return_sim_name()

    if run_sim:

        extPotSim.run_ext_sim(cell_models_folder, measure_coords, 20)

    if plot_sim:
        cell_vmem = np.load(
            join(extPotSim.save_folder, sim_name + '_vmem.npy'))
        cell_tvec = np.load(
            join(extPotSim.save_folder, sim_name + '_tvec.npy'))
        z_rot = [np.pi, 0]

        plotSim = PlotSimulations(
            cellsim_Hallermann_params, monophasic_pulse_params, cell_vmem, cell_tvec)

        plotSim.plot_double_morphology(
            cell_models_folder, z_rot, measure_coords)

    if z == 0 and plot_sim:
        z_rot_len = len(extPotSim.cell_name) + len(f'{extPotSim.x_shift}') + 16
        cell_vmem_rot = np.load(join(extPotSim.save_folder, sim_name[0:z_rot_len] +
                                     f'{np.pi}' + sim_name[z_rot_len + 1:] + '_vmem.npy'))
        cell_tvec_rot = np.load(join(extPotSim.save_folder, sim_name[0:z_rot_len] +
                                     f'{np.pi}' + sim_name[z_rot_len + 1:] + '_tvec.npy'))
        plotSim.plot_double_mem_pot(cell_vmem_rot, cell_tvec_rot)


start = time.time()

task_idx = -1

for z in cell_rot:
    for I in current_amps:
        for pos in elec_positions:
            task_idx += 1
            if not divmod(task_idx, SIZE)[1] == RANK:
                continue

            run_double_hallermann(
                cell_models_folder, measure_coordinates, I, pos, z, True, True)
            print("RANK %d doing task %d" % (RANK, task_idx))

end = time.time()
print(f'Time to execute {end - start} seconds')
