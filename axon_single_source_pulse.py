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

"""

root_folder = os.path.abspath(join(os.path.dirname(__file__), '..'))
cell_models_folder = join(os.path.dirname(__file__), "cell_models")

# print(cell_models_folder)

# model_path = join(cell_models_folder, 'unmyelinated_axon.hoc')

# monophasic_pulse_params["magnitudes"] = [+.001]


# axonsim = ExternalPotentialSim(cellsim_bisc_stick_params)
#
# cell_parameters = axonsim.return_cell(cell_models_folder)
#
# cell = LFPy.Cell(**cell_parameters)
#
# neuron.h('forall insert hh')
#
# axonsim.extra_cellular_stimuli(cell, monophasic_pulse_params)
#
# cell.simulate(rec_vmem=True)
# cell.set_rotation(x=4.729, y=-3.166)
#
# axonsim.plot_cellsim()

cellsim_Hallermann_params['cell_dist_to_top'] = 900

hmsim = ExternalPotentialSim(cellsim_Hallermann_params)

cell_parameters = hmsim.return_cell(cell_models_folder)

neuron.h('forall insert hh')

hmsim.extra_cellular_stimuli(monophasic_pulse_params)

hmsim.plot_cellsim_alt(np.array([0, 83, 300]))


def run_ext_sim(cell, sim_obj, cellsimParams):

    for I in current_amps:

        monophasic_pulse_params['pulse_amp'] = I

        for pos in positions:

            monophasic_pulse_params['positions'] = pos
