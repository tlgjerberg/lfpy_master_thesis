from main import ExternalPotentialSim
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
positions = [np.array([[0, 0, -50], ], dtype=float)]
current_amps = [-1e4]  # uA

axon_measure_idxs = np.array([0, 10, 20, 30, 48])

extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)

extPotSim.run_current_sim(cell_models_folder, monophasic_pulse_params, current_amps,
                          positions, axon_measure_idxs, 10e9, passive=True)
