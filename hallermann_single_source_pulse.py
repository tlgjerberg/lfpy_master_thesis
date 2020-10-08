from main import ExternalPotentialSim
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

extPotSim = ExternalPotentialSim(cellsim_Hallermann_params)


# Test parameters
# current_amps = [1e4, -1e4, -5e3, -4e3, -3e3,  -1e3]  # uA
# positions = [np.array([[200, 0, 700], ], dtype=float),
#              np.array([[201, 130, 699], ], dtype=float),
#              np.array([[-215, 0, 860], ], dtype=float),
#              np.array([[-150, 0, -180], ], dtype=float),
#              np.array([[-20, 0, 400], ], dtype=float)]

current_amps = [-1e4]
positions = [np.array([[200, 0, 700], ], dtype=float)]

coordinates = np.array([[0, 0, -200], [201, 131, 602], [-242, 43, 929]])

# monophasic_pulse_params['stop_time'] = 20
extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params, current_amps,
                      positions, np.array([0, 123, 470]), 20)
