from main import ExternalPotential
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

model_path = join(cell_models_folder, 'unmyelinated_axon.hoc')


def cellParameters(cellsim_params):

    cell_parameters = {
        'morphology': model_path,
        'nsegs_method': "lambda_f",
        'lambda_f': 1000.,
        'v_init': -65,
        'passive': False,
        'dt': cellsim_params['dt'],  # [ms] Should be a power of 2
        'tstart': -cellsim_params['cut_off'],  # [ms] Simulation start time
        'tstop': cellsim_params['tstop'],  # [ms] Simulation end time
        "pt3d": True,
        "extracellular": True,
    }
    return cell_parameters


cell_parameters = cellParameters(cellsim_bisc_stick_params)

cell = LFPy.Cell(**cell_parameters)

simulation = ExternalPotential(cellsim_bisc_stick_params)

simulation.extra_cellular_stimuli(cell, monophasic_pulse_params)
