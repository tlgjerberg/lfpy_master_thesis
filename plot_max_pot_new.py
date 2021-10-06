from neurosim import NeuronSimulation
from plotting import PlotSimulation
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from os.path import join
import sys
import time
from matplotlib.cbook import flatten


cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field/'
cellsim_Hallermann_params['z_rot'] = np.pi
save_folder = 'data/Hallermann_ext_stim/no_field/'
current_amps = [-10.0, -10.5, -11.0, -11.5, -12.0, -12.5, -13.0, -13.5, -14.0, -14.5, -
                15.0, -15.5, -16.0, -16.5, -17.0, -17.5, -18.0, -18.5, -19.0, -19.5, -20.0]


measure_coords = np.array(
    [[0, 0, 0], [-382, 85, 1100], [126, 89, 444], [669, 61, -345]])


segments = ['Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=-100',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=-382',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=226',
            'Hallermann_x_shift=0_z_shift=0_z_rot=3.14_y_rot=0.00_elec_pos=669']

neurosim = NeuronSimulation(cellsim_Hallermann_params)

idx = [0, 471, 104, 125]


def plot_max_potentials(segments, idx):
    v_max_list = []
    for i in range(len(segments)):

        vml = neurosim.consolidate_v_max(segments[i], idx[i])
        # print(vml)
        v_max_list.append(vml)

    # print(v_max_list)
    cell = neurosim.return_cell(cell_models_folder)
    neurosim.create_measure_points(cell, measure_coords)
    plotSim = PlotSimulation(save_folder)
    xlim = [-400, 700]
    ylim = [-400, 1200]
    plotSim.plot_idxs(neurosim.measure_pnts)

    # Plotting
    fig = plt.figure(figsize=[12, 6])
    plotSim.plot_morphology(cell, fig, xlim, ylim, [0.05, 0.1, 0.3, 0.9])
    ax = fig.add_axes([0.45, 0.225, 0.5, 0.6])
    color_list = ['b', 'cyan', 'orange', 'green']
    for c, vm in enumerate(v_max_list):
        ax.plot(current_amps, vm, 'o-', c=color_list[c])

    plt.xlabel('Electode Current Amplitude ($\mu A$)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend(['Soma', 'Apical Dendrite', 'Axon',
                'Axon Terminal', 'Axon Terminal2'])

    # plt.show()
    fig.savefig(join(save_folder,
                     f'v_max_current_dist=100_z_rot={np.pi:.2f}.png'))


plot_max_potentials(segments, idx)
