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
from mpi4py import MPI

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_model'

x_shift = -30
z_shift = 40

measure_coords = np.array(
    [[0, 0, 0], [-382, 85, 1100], [126, 89, 444], [127, 110, 859]])

elec_pos = set_electrode_pos(measure_coords, x_shift, z_shift)


for pos in elec_pos:

    monophasic_pulse_params['positions'] = pos

    plotSim = PlotSimulations(
        cellsim_Hallermann_params, monophasic_pulse_params)

    cell = plotSim.return_cell(cell_models_folder)

    plotSim.create_measure_points(cell, measure_coords)

    plotSim.morphology_3D(cell)
