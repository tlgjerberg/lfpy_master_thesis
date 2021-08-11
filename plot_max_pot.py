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
    """
    Convert compartment index to segment name
    """
    idx_to_sec = {
        '0': 'soma',
        '471': 'apic',
        '104': 'axon',
        '125': 'axon_term',
        '42': 'axon_term2'

    }

    return idx_to_sec.get(num_key, "no index")


x_shift = -40
z_shift = 30

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
# cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim'
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field/'


# current_amps = [1e4, 9e3,  8e3,  7e3, 6e3, 5e3, 4.5e3, 4e3]
# current_amps = [-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4.5e3, -4e3, -3.5e3, -3e3]
current_amps = [-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4.5e3, -4e3]
measure_coords = np.array(
    [[0, 0, 0], [-382, 85, 1100], [126, 89, 444], [127, 110, 859], [669, 61, -345]])


elec_pos = set_electrode_pos(measure_coords, x_shift, z_shift)


measure_keys = ['soma', 'apic', 'axon', 'axon_term', 'axon_term2']
value = []

v_max_sorted = {key: list(value) for key in measure_keys}

z_rot = np.pi
task_idx = -1
for idx, pos in enumerate(elec_pos):
    for I in current_amps:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue

        monophasic_pulse_params['positions'] = pos
        monophasic_pulse_params['pulse_amp'] = I
        cellsim_Hallermann_params['z_rot'] = z_rot

        extPotSim = ExternalPotentialSim(
            cellsim_Hallermann_params, monophasic_pulse_params)

        cell = extPotSim.return_cell(cell_models_folder)

        # Finding the measurement point of index
        extPotSim.create_measure_points(cell, measure_coords[idx])

        elec_idx = str(extPotSim.measure_pnts[0])
        sec = idx_to_sec_conversion(elec_idx)

        extPotSim.extracellular_stimuli(cell)
        extPotSim.run_cell_simulation(cell)
        v_max = extPotSim.max_mem_pot_dict(cell.vmem)
        v_max_sorted[sec].append(v_max[elec_idx])

        cell.__del__()


soma_v_max = COMM.gather(v_max_sorted['soma'], root=0)
apic_v_max = COMM.gather(v_max_sorted['apic'], root=0)
axon_terminal_v_max = COMM.gather(v_max_sorted['axon_term'], root=0)
axon_terminal2_v_max = COMM.gather(v_max_sorted['axon_term2'], root=0)
axon_v_max = COMM.gather(v_max_sorted['axon'], root=0)


if RANK == 0:
    # Sorting and flattening maximum potential into list for each recorded compartment
    soma_v_max_cons = sorted(list(flatten(soma_v_max)), reverse=True)
    apic_v_max_cons = sorted(list(flatten(apic_v_max)), reverse=True)
    axon_terminal_v_max_cons = sorted(
        list(flatten(axon_terminal_v_max)), reverse=True)
    axon_terminal2_v_max_cons = sorted(
        list(flatten(axon_terminal2_v_max)), reverse=True)
    axon_v_max_cons = sorted(list(flatten(axon_v_max)), reverse=True)
    v_max_list = [soma_v_max_cons, apic_v_max_cons, axon_v_max_cons,
                  axon_terminal_v_max_cons, axon_terminal2_v_max_cons]
    current_amps_sorted = sorted(current_amps)

    # Plotting the cell morphology along with recorded compartments
    plotSim = PlotSimulations(
        cellsim_Hallermann_params, monophasic_pulse_params)
    cell = plotSim.return_cell(cell_models_folder)
    plotSim.create_measure_points(cell, measure_coords)
    plotSim.xlim = [-500, 500]
    plotSim.ylim = [-300, 1200]
    plotSim.cell_plot_idxs = plotSim.measure_pnts.astype(
        dtype='int')  # List of measurement points

    plotSim.cell_plot_colors = {idx: [
        'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(plotSim.cell_plot_idxs)}
    fig = plt.figure(figsize=[12, 6])
    plotSim.plot_morphology(cell, fig, [0.05, 0.1, 0.3, 0.9])

    ax = fig.add_axes([0.45, 0.225, 0.5, 0.6])
    color_list = ['b', 'cyan', 'orange', 'green', 'purple']
    for c, vm in enumerate(v_max_list):
        ax.plot(current_amps_sorted, vm, 'o-', c=color_list[c])

    plt.xlabel('Electode Current Amplitude ($\mu A$)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend(['Soma', 'Apical Dendrite', 'Axon',
                'Axon Terminal', 'Axon Terminal2'])
    plt.show()
    fig.savefig(join(extPotSim.save_folder,
                     f'v_max_current_x_shift={x_shift}_z_shift={z_shift}_z_rot={z_rot}.png'))
