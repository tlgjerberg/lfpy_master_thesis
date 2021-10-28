from main import ExternalPotentialSimulation
from plotting import PlotSimulation
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


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_bisc_stick_params['save_folder_name'] = 'data/axon_bisc_field_angle'


def run_axon_angle(cell_models_folder, measure_coords, I, pos, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    monophasic_pulse_params['stop_time'] = 10.0
    cellsim_bisc_stick_params['cell_dist_to_top'] = -500

    extPotSim = ExternalPotentialSimulation(
        cellsim_bisc_stick_params, monophasic_pulse_params)
    cell = extPotSim.return_cell(cell_models_folder)

    dV_pot_dict = None

    if run_sim:

        elec_positions = set_electrode_pos(measure_coords)
        extPotSim.run_ext_sim(cell, measure_coords)
        dV_pot_dict = extPotSim.dV_pot_dict(cell.vmem)

        return dV_pot_dict

    else:
        print('No simulation run!')

    if plot_sim:

        extPotSim.import_data()
        extPotSim.create_measure_points(cell, measure_coords)
        extPotSim.print_measure_points(cell)
        plotSim = PlotSimulation(extPotSim.save_folder)
        fig = plt.figure()
        plotSim.plot_idxs(extPotSim.measure_pnts)
        plotSim.legend_list = ['Compartment 24',
                               'Compartment 0', 'Compartment 48']
        plotSim.plot_membrane_potential(
            fig, extPotSim.cell_tvec, extPotSim.cell_vmem, extPotSim.tstop)

        sim_name = extPotSim.return_sim_name()
        fig.savefig(join(extPotSim.save_folder,
                         f'mem_pot_{sim_name}.png'), dpi=300)
        fig_morph = plt.figure()
        plotSim.plot_morphology(cell, fig_morph, [-100, 2000], [-2000, 2000])
        plotSim.draw_electrode(extPotSim.x0, extPotSim.y0,
                               extPotSim.z0, extPotSim.electrode_radii)
        plt.show()
    else:
        print('No plots generated!')

    # return v_max


measure_coords = np.array([[0, 0, -500], [0, 0, 0], [0, 0, 500]])


elec_positions = np.array([[2000, 0, 0],
                           [2000 * np.cos(np.pi / 4), 0, 2000 *
                            np.sin(np.pi / 4)],
                           [2000 * np.cos(np.pi / 3), 0, 2000 *
                            np.sin(np.pi / 3)],
                           [2000 * np.cos(np.pi / 6), 0, 2000 *
                            np.sin(np.pi / 6)],
                           [2000 * np.cos(np.pi / 2), 0, 2000 *
                            np.sin(np.pi / 2)],
                           [2000 * np.cos(3 * np.pi / 2), 0,
                            2000 * np.sin(3 * np.pi / 2)],
                           [2000 * np.cos(7 * np.pi / 4), 0,
                            2000 * np.sin(7 * np.pi / 4)],
                           [0, 0, 2000]], dtype=float)

# elec_positions = np.array([[2000, 0, 500],
#                            [2000 * np.cos(np.pi / 4), 0, 2500 *
#                             np.sin(np.pi / 4)],
#                            [2000 * np.cos(np.pi / 2), 0, 2500 *
#                             np.sin(np.pi / 2)],
#                            [2000 * np.cos(np.pi / 6), 0, 2500 *
#                             np.sin(np.pi / 6)],
#                            [2000 * np.cos(np.pi / 3), 0, 2500 *
#                             np.sin(np.pi / 3)]], dtype=float)

dV_pot_dict = None
I = -5e4  # uA

task_idx = -1
for pos in elec_positions:

    task_idx += 1
    if not divmod(task_idx, SIZE)[1] == RANK:
        continue

    dV_pot_dict = run_axon_angle(
        cell_models_folder, measure_coords, I, pos, True)
    print('Electrode position: ', pos)
    print('maximum potential: ', dV_pot_dict)

    print("RANK %d doing task %d" % (RANK, task_idx))


# if v_max is not None:
#     print(f'RANK {RANK} has reached gather')
#     v_max_top = COMM.gather(v_max['0'], root=0)
#     v_max_mid = COMM.gather(v_max['24'], root=0)
#     v_max_bot = COMM.gather(v_max['48'], root=0)
#
# if RANK == 0:
#     print(v_max_top)
#     v_max_top_cons = sorted(list(flatten(v_max_top)), reverse=True)
#     v_max_mid_cons = sorted(list(flatten(v_max_mid)), reverse=True)
#     v_max_bot_cons = sorted(list(flatten(v_max_bot)), reverse=True)
#     v_max_list = [v_max_top_cons, v_max_mid_cons, v_max_bot_cons]
#     radians = [0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2]
#     radians_sorted = sorted(radians)
