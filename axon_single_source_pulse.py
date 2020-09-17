from main import ExternalPotentialSim
from parameters import (monophasic_pulse_params,
                        cellsim_bisc_stick_params, cellsim_Hallermann_params)
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


def run_ext_sim(cellsimParams, elec_params, I, positions, passive=False):

    extPotSim = ExternalPotentialSim(cellsimParams)
    extPotSim.return_cell(cell_models_folder)

    # Neuron activation after cell object has been created
    if not passive:
        neuron.h('forall insert hh')

    for I in current_amps:

        monophasic_pulse_params['pulse_amp'] = I

        for pos in positions:

            monophasic_pulse_params['positions'] = pos
            extPotSim.extra_cellular_stimuli(monophasic_pulse_params)
            extPotSim.plot_cellsim(np.array([0, 83, 300]))


# Test parameters
current_amps = [1e4, -1e4, 5e3]  # uA
# positions = [np.array([[200, 0, -40], ], dtype=float),
#              np.array([[200, 0, 0], ], dtype=float),
#              np.array([[-125, 0, -880], ], dtype=float),
#              np.array([[-230, 0, 175], ], dtype=float)]
# current_amps = [-1e4]  # uA
positions = [np.array([[210, 0, -40], ], dtype=float)]
cellsim_Hallermann_params['cell_dist_to_top'] = 900


run_ext_sim(cellsim_Hallermann_params,
            monophasic_pulse_params, current_amps, positions)

# extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)
# extPotSim.plot_axialCurrent()
