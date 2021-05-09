from main import ExternalPotentialSim
from plotting import PlotSimulations
from parameters import (monophasic_pulse_params, cellsim_bisc_stick_params)
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_bisc_stick_params['save_folder_name'] = "data/axon_bisc_current_snap"
positions = np.array([0, 0, -50], dtype=float)

monophasic_pulse_params['pulse_amp'] = -1e4  # uA
monophasic_pulse_params['positions'] = positions

measure_coordinates = np.array(
    [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])
z = np.pi

extPotSim = ExternalPotentialSim(
    cellsim_bisc_stick_params, monophasic_pulse_params)


extPotSim.run_current_sim(
    cell_models_folder, measure_coordinates, passive=True)

cell_vmem = None

cell_tvec = np.load(
    join(extPotSim.save_folder, f'axon_x_shift=0_z_rot=0.00_{-10000.0}mA_elec_pos={0.0}_{0.0}_{-50.0}_tvec' + '.npy'))
cell_imem = np.load(
    join(extPotSim.save_folder, f'axon_x_shift=0_z_rot=0.00_{-10000.0}mA_elec_pos={0.0}_{0.0}_{-50.0}_imem' + '.npy'))

plotSim = PlotSimulations(cellsim_bisc_stick_params,
                          monophasic_pulse_params, cell_vmem, cell_tvec, cell_imem)
cell = plotSim.return_cell(cell_models_folder)
plotSim.plot_currents(cell, measure_coordinates, [0.05, 0.05, 0.3, 0.90])
