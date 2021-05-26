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
chiral_morphology
"""

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_bisc_stick_params['save_folder_name'] = 'data/axon_bisc_field_angle'


def run_axon_angle(cell_models_folder, measure_coords, I, pos, y_rot, x_shift, z_shift, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_bisc_stick_params['y_rot'] = y_rot
    cellsim_bisc_stick_params['x_shift'] = x_shift
    cellsim_bisc_stick_params['cell_dist_to_top'] = z_shift

    extPotSim = ExternalPotentialSim(
        cellsim_bisc_stick_params, monophasic_pulse_params)

    if run_sim:
        elec_positions = set_electrode_pos(measure_coords)

        extPotSim.run_ext_sim(cell_models_folder, measure_coords, z)

    else:
        print('No simulation run!')

    if plot_sim:

        cell_tvec = np.load(
            join(extPotSim.save_folder, f'axon_x_shift={x_shift}_z_shift={z_shift}_z_rot=0.00_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_tvec' + '.npy'), allow_pickle=True)
        cell_vmem = np.load(
            join(extPotSim.save_folder, f'axon_x_shift={x_shift}_z_shift={z_shift}_z_rot=0.00_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_vmem' + '.npy'), allow_pickle=True)
        plotSim = PlotSimulations(
            cellsim_bisc_stick_params, monophasic_pulse_params, cell_vmem, cell_tvec)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim(measure_coords, z, [
                             0.05, 0.05, 0.3, 0.90], [-1200, 1200], [-1200, 100], field=True)

    else:
        print('No plots generated!')


measure_coords = np.array(
    [[0, 0, -600], [0, 0, -1100], [0, 0, 0], [-600, 0, -1100], [600, 0, -1100], [424, 0, -424]])
elec_positions = np.array([[0, 0, 0]], dtype=float)
pos = elec_positions[0]
z = np.pi
task_idx = -1
I = -1e4  # uA
cell_rot = [0, np.pi / 2, ]
x_shift = [0, -500]
z_shift = [-1100, -1600]
for y_rot in cell_rot:

    for x in x_shift:

        for z in z_shift:
            task_idx += 1
            if not divmod(task_idx, SIZE)[1] == RANK:
                continue

            run_axon_angle(cell_models_folder, measure_coords,
                           I, pos, y_rot, x, z, True, True)

            print("RANK %d doing task %d" % (RANK, task_idx))
