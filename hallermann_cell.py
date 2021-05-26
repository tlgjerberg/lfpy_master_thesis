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

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_ext_stim/no_field'


cell_rot = [np.pi, (1 / 2) * np.pi]
for z_rot in cell_rot:
    cellsim_Hallermann_params['z_rot'] = z_rot
    plotSim = PlotSimulations(
        cellsim_Hallermann_params, monophasic_pulse_params)
    cell = plotSim.return_cell(cell_models_folder)
    # plotSim.morphology_3D(cell)
    fig = plt.figure()
    plotSim.xlim = [-500, 500]
    plotSim.ylim = [-300, 1200]
    plotSim.plot_morphology(cell, fig)
