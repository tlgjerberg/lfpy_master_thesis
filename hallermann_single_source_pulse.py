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
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim'
# cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'

# current_amps = [1e4, -1e4, -7e3, 7e3]
current_amps = [1e4, -1e4]


measure_coordinates = np.array(
    [[0, 0, 0], [127, 126, 866], [-393, 39, 1101], [123, 100, 443]])

elec_positions = set_electrode_pos(measure_coordinates)


def run_hallermann(cell_models_folder, measure_coords, I, pos, z=np.pi, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_Hallermann_params['z_rot'] = z

    extPotSim = ExternalPotentialSim(
        cellsim_Hallermann_params, monophasic_pulse_params)

    if run_sim:

        extPotSim.run_ext_sim(cell_models_folder, measure_coords, 20, z)

    if plot_sim:
        cell_vmem = np.load(join(extPotSim.save_folder,
                                 f'Hallermann_x_shift=0_z_rot={z:.2f}_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_vmem.npy'))
        cell_tvec = np.load(join(extPotSim.save_folder,
                                 f'Hallermann_x_shift=0_z_rot={z:.2f}_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_tvec.npy'))
        plotSim = PlotSimulations(
            cellsim_Hallermann_params, monophasic_pulse_params, cell_vmem, cell_tvec)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(cell, measure_coords, z, [0.05, 0.05, 0.3, 0.90])


start = time.time()
z = np.pi
task_idx = -1
for I in current_amps:
    for pos in elec_positions:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue

        run_hallermann(cell_models_folder,
                       measure_coordinates, I, pos, z, True, True)
        print("RANK %d doing task %d" % (RANK, task_idx))


end = time.time()
print(f'Time to execute {end - start} seconds')
