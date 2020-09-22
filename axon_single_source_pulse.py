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

"""

Improve plotting

Add method for running over a set of electrode positions and current amplitudes

Create plot of Hallermann pyrmidal cell (4 figs?) showing electrode near soma,
dendrite and axon to demonstrate if activation is possible and at what current
amplitude.

Stimulate using current of 10 muA and go lower after (Histed et Al)

Create figure of axial current over distance along axon. Use simple stick model.

Improve plot_cellsim_alt for easy reading and page formatting.


"""

cell_models_folder = join(os.path.dirname(__file__), "cell_models")

# Test parameters
# current_amps = [1e4, -1e4, 5e3]  # uA
positions = [np.array([[200, 0, -40], ], dtype=float),
             np.array([[200, 0, 0], ], dtype=float),
             np.array([[-125, 0, -880], ], dtype=float),
             np.array([[-230, 0, 175], ], dtype=float)]
current_amps = [-1e4]  # uA
# positions = [np.array([[210, 0, 700], ], dtype=float)]

axon_measure_idxs = np.array([0, 20, 48])

extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)

extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params, current_amps,
                      positions, axon_measure_idxs, 200)
