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


def run_axon_angle(cell_models_folder, measure_coords, I, pos, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_bisc_stick_params['cell_dist_to_top'] = -500
    cellsim_bisc_stick_params['tstop'] = 20
    print('RANK', RANK, 'electrode position', pos)

    extPotSim = ExternalPotentialSim(
        cellsim_bisc_stick_params, monophasic_pulse_params)

    if run_sim:
        elec_positions = set_electrode_pos(measure_coords)

        v_max = extPotSim.run_ext_sim(cell_models_folder, measure_coords)

    else:
        print('No simulation run!')

    if plot_sim:

        cell_tvec = np.load(
            join(extPotSim.save_folder, extPotSim.sim_name + '_tvec.npy'))
        cell_vmem = np.load(
            join(extPotSim.save_folder, extPotSim.sim_name + '_vmem.npy'))

        plotSim = PlotSimulations(
            cellsim_bisc_stick_params, monophasic_pulse_params, cell_vmem, cell_tvec)
        cell = plotSim.return_cell(cell_models_folder)
        plotSim.plot_cellsim_angle(
            measure_coords, [0, 0, 1, 1], field=True)
    else:
        print('No plots generated!')

    return v_max


measure_coords = np.array([[0, 0, -10], [0, 0, -500], [0, 0, 500]])

# elec_positions = np.array([[600, 0, -555],
#                            [600 * np.cos(np.pi / 4), 0, -555 +
#                             600 * np.sin(np.pi / 4)],
#                            [600 * np.cos(np.pi / 3), 0, -555 +
#                             600 * np.sin(np.pi / 3)],
#                            [600 * np.cos(np.pi / 6), 0, -555 +
#                             600 * np.sin(np.pi / 6)],
#                            [600 * np.cos(np.pi / 2), 0, -555 +
#                             600 * np.sin(np.pi / 2)],
#                            [0, 0, 600]], dtype=float)
# elec_positions = np.array([[1000, 0, -555],
#                            [1000 * np.cos(np.pi / 4), 0, -555 +
#                             1000 * np.sin(np.pi / 4)],
#                            [1000 * np.cos(np.pi / 3), 0, -555 +
#                             1000 * np.sin(np.pi / 3)],
#                            [1000 * np.cos(np.pi / 6), 0, -555 +
#                             1000 * np.sin(np.pi / 6)],
#                            [1000 * np.cos(np.pi / 2), 0, -555 +
#                             1000 * np.sin(np.pi / 2)],
#                            [0, 0, 600]], dtype=float)
elec_positions = np.array([[2000, 0, -10],
                           [2000 * np.cos(np.pi / 4), 0, -10 +
                            2000 * np.sin(np.pi / 4)],
                           [2000 * np.cos(np.pi / 3), 0, -10 +
                            2000 * np.sin(np.pi / 3)],
                           [2000 * np.cos(np.pi / 6), 0, -10 +
                            2000 * np.sin(np.pi / 6)],
                           [2000 * np.cos(np.pi / 2), 0, -10 +
                            2000 * np.sin(np.pi / 2)],
                           [2000 * np.cos(3 * np.pi / 2), 0, -10 +
                            2000 * np.sin(3 * np.pi / 2)],
                           [2000 * np.cos(7 * np.pi / 4), 0, -10 +
                            2000 * np.sin(7 * np.pi / 4)]], dtype=float)

measure_keys = ['0', '25', '48']
value = []

# v_max_sorted = {key: list(value) for key in measure_keys}
# print(v_max_sorted)
# v_max = None
I = -1e4  # uA

task_idx = -1
for pos in elec_positions:

    task_idx += 1
    if not divmod(task_idx, SIZE)[1] == RANK:
        continue

    v_max = run_axon_angle(
        cell_models_folder, measure_coords, I, pos, True, True)

    print("RANK %d doing task %d" % (RANK, task_idx))

# print('v_max', v_max)
# print(v_max_sorted)
# if v_max is not None:
#     print(f'RANK {RANK} has reached gather')
#     v_max_top = COMM.gather(v_max['0'], root=0)
#     v_max_mid = COMM.gather(v_max['25'], root=0)
#     v_max_bot = COMM.gather(v_max['48'], root=0)
# #
# if RANK == 0:
#     print(v_max_top)
# #     v_max_top_cons = sorted(list(flatten(v_max_top)), reverse=True)
#     v_max_mid_cons = sorted(list(flatten(v_max_mid)), reverse=True)
#     v_max_bot_cons = sorted(list(flatten(v_max_bot)), reverse=True)
#     v_max_list = [v_max_top_cons, v_max_mid_cons, v_max_bot_cons]
#     radians = [0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2]
#     radians_sorted = sorted(radians)
#
#     plotSim = PlotSimulations(
#         cellsim_Hallermann_params, monophasic_pulse_params)
#     cell = plotSim.return_cell(cell_models_folder)
#     plotSim.create_measure_points(cell, measure_coords)
#     plotSim.xlim = [-1050, 1050]
#     plotSim.ylim = [-1200, 800]
#     plotSim.cell_plot_idxs = plotSim.measure_pnts.astype(
#         dtype='int')  # List of measurement points
#
#     plotSim.cell_plot_colors = {idx: [
#         'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(plotSim.cell_plot_idxs)}
#     fig = plt.figure(figsize=[12, 6])
#     plotSim.plot_morphology(cell, fig, [0.05, 0.1, 0.3, 0.9])
#
#     ax = fig.add_axes([0.45, 0.225, 0.5, 0.6])
#     color_list = ['b', 'cyan', 'orange', 'green', 'purple']
#     for c, vm in enumerate(v_max_list):
#         ax.plot(current_amps_sorted, vm, 'o-', c=color_list[c])
#
#     plt.xlabel('Electode Current Amplitude ($\mu A$)')
#     plt.ylabel('Membrane Potential (mV)')
#     plt.legend(['Soma', 'Apical Dendrite', 'Axon Terminal', 'Axon'])
#     plt.show()
#     fig.savefig(join(extPotSim.save_folder,
#                      f'v_max_current_x_shift={x_shift}_z_shift={z_shift}_z_rot={z_rot}.png'))
