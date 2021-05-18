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
from matplotlib.cbook import flatten

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


def idx_to_sec_conversion(num_key='0'):

    idx_to_sec = {
        '0': 'soma',
        '392': 'apic',
        '352': 'axon',
        '521': 'axon_term'
    }

    return idx_to_sec.get(num_key, "no index")


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
# cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim'
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'


# current_amps = [1e4, 9e3,  8e3,  7e3, 6e3, 5e3, 4.5e3, 4e3]
# current_amps = [-1e4, -9e3, -8e3, -7e3]
current_amps = [-6e3, -5e3, -4.5e3, -4e3]

measure_coords = np.array(
    [[0, 0, 0], [-393, 80, 1101], [123, 90, 443], [127, 126, 866]])

elec_pos = set_electrode_pos(measure_coords, -30, 40)


measure_keys = ['soma', 'apic', 'axon', 'axon_term']
value = []
# v_max_sorted = dict.fromkeys(measure_keys, [])

v_max_sorted = {key: list(value) for key in measure_keys}

task_idx = -1
for idx, pos in enumerate(elec_pos):
    for I in current_amps:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue

        monophasic_pulse_params['positions'] = pos
        monophasic_pulse_params['pulse_amp'] = I
        extPotSim = ExternalPotentialSim(
            cellsim_Hallermann_params, monophasic_pulse_params)

        cell = extPotSim.return_cell(cell_models_folder)
        extPotSim.create_measure_points(cell, measure_coords[idx])

        if pos[1] == 0:
            pos = measure_coords[0]

        elec_idx = str(cell.get_closest_idx(pos[0], pos[1], pos[2]))
        print('RANK', RANK, 'elec_idx', elec_idx)
        sec = idx_to_sec_conversion(elec_idx)
        print('RANK', RANK, 'sec', sec)

        # extPotSim.print_measure_points(cell)

        extPotSim.extracellular_stimuli(cell)
        extPotSim.run_cell_simulation(cell)
        v_max = extPotSim.find_max_mem_pot(cell.vmem)
        print(v_max_sorted[sec])
        print('v_max', v_max)
        print('v_max_elec', v_max[elec_idx])
        v_max_sorted[sec].append(v_max[elec_idx])
        print('v_max_sorted', v_max_sorted)
        cell.__del__()
        print("RANK %d doing task %d" % (RANK, task_idx))


# if RANK == 0:
#     print('RANK: ', RANK)

soma_v_max = COMM.gather(v_max_sorted['soma'], root=0)
apic_v_max = COMM.gather(v_max_sorted['apic'], root=0)
axon_terminal_v_max = COMM.gather(v_max_sorted['axon_term'], root=0)
axon_v_max = COMM.gather(v_max_sorted['axon'], root=0)

# else:
#     print('RANK: ', RANK)
#     if sec == 'soma':
#         soma_v_max = COMM.gather(v_max_sorted['soma'], root=0)
#     elif sec == 'apic':
#         apic_v_max = COMM.gather(v_max_sorted['apic'], root=0)
#     elif sec == 'axon':
#         axon_terminal_v_max = COMM.gather(v_max_sorted['axon_term'], root=0)
#     elif sec == 'axon_term':
#         axon_v_max = COMM.gather(v_max_sorted['axon'], root=0)


if RANK == 0:
    soma_v_max_cons = sorted(list(flatten(soma_v_max)))
    apic_v_max_cons = sorted(list(flatten(apic_v_max)))
    axon_terminal_v_max_cons = sorted(list(flatten(axon_terminal_v_max)))
    axon_v_max_cons = sorted(list(flatten(axon_v_max)))
    current_amps_sorted = sorted(current_amps)

    plt.plot(current_amps_sorted, soma_v_max_cons)
    plt.plot(current_amps_sorted, axon_terminal_v_max_cons)
    plt.plot(current_amps_sorted, apic_v_max_cons)
    plt.plot(current_amps_sorted, axon_v_max_cons)
    plt.xlabel('Electode Current Amplitude ($\mu A$)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend(['Soma', 'Axon Terminal', 'Apical Dendrite', 'Axon'])
    plt.show()
    plt.savefig(
        join(extPotSim.save_folder, f'v_max_current.png'), dpi=300)
