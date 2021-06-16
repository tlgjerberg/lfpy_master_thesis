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

x_shift = -100
z_shift = 0
# x_shift = -30
# z_shift = 40
# x_shift = -40
# z_shift = 30
# x_shift = 100
# z_shift = 0
# x_shift = 0
# z_shift = 100


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'


# current_amps = [1e4, 9e3,  8e3,  7e3, 6e3, 5e3, 4.5e3, 4e3]
# current_amps = [-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4.5e3, -4e3]
# current_amps = [-3.5e3, -3e3]
# current_amps = [-2e4, -1.9e4, -1.8e4, -1.7e4, -
#                 1.6e4, -1.5e4, -1.4e4, -1.3e4, -1.2e4, -1.1e4, -1e4]
current_amps = [-2e4, -1.9e4, -1.8e4, -1.7e4, -1.65e4, -1.6e4, -1.5e4]
# current_amps = [-1.4e4, -1.3e4, -1.2e4, -1.1e4, -1e4]


measure_coords = np.array(
    [[0, 0, 0], [-382, 85, 1100], [126, 89, 444], [669, 61, -345]])


elec_pos = set_electrode_pos(measure_coords, x_shift, z_shift)


def run_hallermann(cell_models_folder, measure_coords, I, pos, z=np.pi, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_Hallermann_params['z_rot'] = z
    # monophasic_pulse_params['stop_time'] = 2.2

    extPotSim = ExternalPotentialSim(
        cellsim_Hallermann_params, monophasic_pulse_params)

    if run_sim:

        extPotSim.run_ext_sim(cell_models_folder, measure_coords)

    if plot_sim:
        extPotSim.return_sim_name()
        cell_vmem = np.load(
            join(extPotSim.save_folder, extPotSim.sim_name + '_vmem.npy'))
        cell_tvec = np.load(
            join(extPotSim.save_folder, extPotSim.sim_name + '_tvec.npy'))
        plotSim = PlotSimulations(
            cellsim_Hallermann_params, monophasic_pulse_params, cell_vmem, cell_tvec)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(measure_coords, [0.05, 0.05, 0.3, 0.90])


start = time.time()
z = np.pi
task_idx = -1
for I in current_amps:
    for pos in elec_pos:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue
        run_hallermann(
            cell_models_folder, measure_coords, I, pos, z, True, True)
        print("RANK %d doing task %d" % (RANK, task_idx))


end = time.time()
print(f'Time to execute {end - start} seconds')
