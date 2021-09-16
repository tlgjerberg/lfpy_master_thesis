from neurosim import NeuronSimulation
from main import ExternalPotentialSimulation
from plotting import PlotSimulation
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
from comp_idx_to_sec import idx_to_sec_conversion
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from os.path import join
import sys
import time
from mpi4py import MPI
from matplotlib.cbook import flatten
from glob import glob


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field/'
cellsim_Hallermann_params['z_rot'] = np.pi
save_folder = 'data/Hallermann_ext_stim/no_field/'
current_amps = [-1, -1.1, -1.2, -1.3, -1.4, -
                1.5, -1.6, -1.65, -1.7, -1.8, -1.9, -2]


measure_coords = np.array(
    [[0, 0, 0], [-382, 85, 1100], [126, 89, 444], [127, 110, 859], [669, 61, -345]])


def consolidate_v_max(sim_name, idx):
    """
    Loads a set of memembrane potential simulations and consolidates a list of
    all maximum potentials.

    Input:
    segement: String giving the path to a set of simulation.
    idx: index of the compartment of interest.

    Returns:
    vml: A list of maximum potentials.
    """

    vmem_list = sorted(glob(
        join('data/Hallermann_ext_stim/no_field', sim_name + '*_vmem.npy')))
    vml = []
    for vmem in vmem_list:
        cell_vmem = np.load(vmem)
        v_max = np.max(cell_vmem[idx])
        vml.append(v_max)
    return vml


segments = ['Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=-100',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=-482',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=226',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=669']

idx = [0, 471, 104, 125, 42]

v_max_list = []
for i in range(len(segments)):

    vml = consolidate_v_max(segments[i], idx[i])
    v_max_list.append(vml)


neurosim = NeuronSimulation(cellsim_Hallermann_params)
cell = neurosim.return_cell(cell_models_folder)
neurosim.create_measure_points(cell, measure_coords)
plotSim = PlotSimulation(save_folder)
xlim = [-500, 500]
ylim = [-300, 1200]
plotSim.plot_idxs(neurosim.measure_pnts)

# Plotting
fig = plt.figure(figsize=[12, 6])
plotSim.plot_morphology(cell, fig, xlim, ylim, [0.05, 0.1, 0.3, 0.9])
ax = fig.add_axes([0.45, 0.225, 0.5, 0.6])
color_list = ['b', 'cyan', 'orange', 'purple', 'green']
for c, vm in enumerate(v_max_list):
    ax.plot(current_amps, vm, 'o-', c=color_list[c])

plt.xlabel('Electode Current Amplitude ($\mu A$)')
plt.ylabel('Membrane Potential (mV)')
plt.legend(['Soma', 'Apical Dendrite', 'Axon',
            'Axon Terminal', 'Axon Terminal2'])

# plt.show()
fig.savefig(join(save_folder,
                 f'v_max_current_dist=100_z_rot={np.pi:.2f}.png'))
