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

root_folder = os.path.abspath(join(os.path.dirname(__file__), '..'))
cell_models_folder = join(os.path.dirname(__file__), "cell_models")

print(cell_models_folder)

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

hmsim = ExternalPotentialSim(cellsim_Hallermann_params)

cell_parameters = hmsim.return_cell(cell_models_folder)

cell = LFPy.Cell(**cell_parameters)

neuron.h('forall insert hh')

hmsim.extra_cellular_stimuli(cell, monophasic_pulse_params)

cell.simulate(rec_vmem=True)
cell.set_rotation(x=4.729, y=-3.166, z=-3)

hmsim.plot_cellsim()
