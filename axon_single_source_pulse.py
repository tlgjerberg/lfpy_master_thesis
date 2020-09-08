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

root_folder = os.path.abspath(join(os.path.dirname(__file__), '..'))
cell_models_folder = join(os.path.dirname(__file__), "cell_models")

# model_path = join(cell_models_folder, 'unmyelinated_axon.hoc')


# def cellParameters(cellsim_params):
#
#     cell_parameters = {
#         'morphology': model_path,
#         'nsegs_method': "lambda_f",
#         'lambda_f': 1000.,
#         'v_init': -65,
#         'passive': False,
#         'dt': cellsim_params['dt'],  # [ms] Should be a power of 2
#         'tstart': -cellsim_params['cut_off'],  # [ms] Simulation start time
#         'tstop': cellsim_params['tstop'],  # [ms] Simulation end time
#         "pt3d": True,
#         "extracellular": True,
#     }
#     return cell_parameters


extsim = ExternalPotentialSim(cellsim_bisc_stick_params)

cell_parameters = extsim.return_cell(cell_models_folder)

cell = LFPy.Cell(**cell_parameters)

extsim.extra_cellular_stimuli(cell, monophasic_pulse_params)

cell.simulate(rec_vmem=True)
# cell.set_rotation(x=4.729, y=-3.166)

extsim.plot_cellsim()
