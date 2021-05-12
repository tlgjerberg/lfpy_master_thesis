from main import ExternalPotentialSim
from plotting import PlotSimulations
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from os.path import join
import sys
import time
from mpi4py import MPI

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
# cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim'
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'


# current_amps = [1e4, 9e3,  8e3,  7e3, 6e3, 5e3, 4.5e3, 4e3]
current_amps = [-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4.5e3, -4e3]

measure_coords = np.array(
    [[0, 0, 0], [-393, 80, 1101], [127, 126, 866], [123, 90, 443]])

elec_pos = set_electrode_pos(measure_coords, -30, 40)


measure_keys = [0, 392, 521, 552]
v_max_sorted = dict.fromkeys(measure_keys, [])

task_idx = -1
for pos in elec_pos:
    for I in current_amps:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue

        monophasic_pulse_params['positions'] = pos
        monophasic_pulse_params['pulse_amp'] = I
        extPotSim = PlotSimulations(
            cellsim_Hallermann_params, monophasic_pulse_params)

        cell = extPotSim.return_cell(cell_models_folder)
        extPotSim.create_measure_points(cell, measure_coords)

        if pos[1] == 0:
            pos = measure_coords[0]

        elec_idx = cell.get_closest_idx(pos[0], pos[1], pos[2])

        # measure_idx = np.extract(extPotSim.measure_pnts ==
        #                          elec_idx, extPotSim.measure_pnts)

        # extPotSim.print_measure_points(cell)
        extPotSim.return_segment_coords(cell)

        extPotSim.extracellular_stimuli(cell)
        extPotSim.run_cell_simulation(cell)
        v_max = extPotSim.find_max_mem_pot(cell.vmem)

        v_max_sorted[elec_idx].append(v_max[str(elec_idx)])

        cell.__del__()

soma_v_max = v_max_sorted['0']
apic_v_max = v_max_sorted['392']
axon_terminal_v_max = v_max_sorted['521']
axon_v_max = v_max_sorted['392']
# start = time.time()
# z = np.pi
# v_max = None
# task_idx = -1
# for I in current_amps:
#     for pos in elec_pos:
#         task_idx += 1
#         if not divmod(task_idx, SIZE)[1] == RANK:
#             continue
#
#         v_max = plot_max_pot(
#             cell_models_folder, measure_coords, I, pos, z, True, True)
#         print("RANK %d doing task %d" % (RANK, task_idx))
#
# v_max = COMM.gather(v_max, root=0)
#
# if RANK == 0:
#     extPotSim = ExternalPotentialSim(
#         cellsim_Hallermann_params, monophasic_pulse_params)
#     v_max = [i for i in v_max if i]
#     v_max = np.array(v_max)
#     np.save(join(extPotSim.save_folder,
#                  f'max_pot_elec_pos={elec_pos[0][0]}_{elec_pos[0][1]}_{elec_pos[0][2]}.npy'), v_max)
#
#
# end = time.time()
# print(f'Time to execute {end - start} seconds')
#
# v_max = np.load(join(extPotSim.save_folder,
#                      f'max_pot_elec_pos={elec_pos[0][0]}_{elec_pos[0][1]}_{elec_pos[0][2]}.npy'))
#
#
# plt.figure()
#
# for ep in elec_pos:
#     print(ep)
#     v_max = np.load(join(extPotSim.save_folder,
#                          f'max_pot_elec_pos={ep[0]}_{ep[1]}_{ep[2]}.npy'))
#
#     plt.plot(current_amps, v_max, 'o-')
#
# plt.xlabel('Electode Current Amplitude ($\mu A$)')
# plt.ylabel('Membrane Potential (mV)')
# plt.legend(['Soma', 'Axon Terminal', 'Apical Dendrite', 'Axon'])
# plt.show()
# plt.savefig(
#     join(extPotSim.save_folder, f'v_max_current_elec_pos={elec_pos[0][0]}_{elec_pos[0][1]}_{elec_pos[0][2]}.png'), dpi=300)
