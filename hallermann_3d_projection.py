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

elec_positions = np.array([[-50, 0, 0],
                           [60, 126, 659],
                           [-301, 39, 879]])


# monophasic_pulse_params['positions'] = [-127, -126, 866]
# monophasic_pulse_params['positions'] = [-127, -126, 866]


measure_coordinates = np.array(
    [[0, 0, 0], [-127, -126, 866], [-393, 39, 1101]])

plotSim = PlotSimulations(cellsim_Hallermann_params, monophasic_pulse_params)

cell = plotSim.return_cell(cell_models_folder)

plotSim.create_measure_points(cell, measure_coordinates)

plotSim.morphology_3D(cell)
