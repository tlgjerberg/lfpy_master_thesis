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


def run_dv_pos(cellsimParams, elec_params, I, positions, measure_idxs, passive=False):
    extPotSim = ExternalPotentialSim(cellsimParams)
    extPotSim.return_cell(cell_models_folder)

    # Neuron activation after cell object has been created
    if not passive:
        neuron.h('forall insert hh')

    elec_abs_dists = np.zeros((len(positions), 3))
    ss_pot = np.zeros(len(positions))
    print(elec_abs_dists)

    for I in current_amps:

        monophasic_pulse_params['pulse_amp'] = I

        for idx, pos in enumerate(positions):

            monophasic_pulse_params['positions'] = pos
            extPotSim.extra_cellular_stimuli(monophasic_pulse_params)
            # extPotSim.plot_cellsim(measure_idxs)
            extPotSim.run_cell_simulation()

            elec_abs_dists[idx], ss_pot[idx] = extPotSim.record_dist_to_electrode(
                measure_idxs)

    extPotSim.plot_potentialVdistance(elec_abs_dists[:, 0], ss_pot)


# Test parameters
# current_amps = [1e4, -1e4, 5e3]  # uA
positions = [np.array([[200, 0, -40], ], dtype=float),
             np.array([[200, 0, 0], ], dtype=float),
             np.array([[-125, 0, -880], ], dtype=float),
             np.array([[-230, 0, 175], ], dtype=float)]
current_amps = [-1e4]  # uA
# positions = [np.array([[210, 0, 700], ], dtype=float)]

# run_ext_sim(cellsim_bisc_stick_params,
#             monophasic_pulse_params, current_amps, positions, np.array([0, 20, 48]))
axon_measure_idxs = np.array([0, 20, 48])

# extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)
# extPotSim.return_cell(cell_models_folder)
# extPotSim.extra_cellular_stimuli(monophasic_pulse_params)
# extPotSim.plot_currentVdistance(axon_measure_idxs)

run_dv_pos(cellsim_bisc_stick_params,
           monophasic_pulse_params, current_amps, positions, np.array([0, 20, 48]))
