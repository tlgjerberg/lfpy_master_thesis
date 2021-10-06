from plotting import PlotSimulation
from neurosim import NeuronSimulation
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
from os.path import join
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
import time

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'
save_folder = cellsim_Hallermann_params['save_folder_name']


cell_rot = [np.pi, (1 / 2) * np.pi]
for z_rot in cell_rot:
    cellsim_Hallermann_params['z_rot'] = z_rot

    neuroSim = NeuronSimulation(cellsim_Hallermann_params)
    cell = neuroSim.return_cell(cell_models_folder)
    plotSim = PlotSimulation(save_folder)

    fig = plt.figure()
    xlim = [-600, 600]
    ylim = [-400, 1200]
    plotSim.plot_morphology(cell, fig, xlim, ylim)
    fig.savefig(join(save_folder, f'hallermann_cell_z_rot={z_rot:.2f}.png'))
