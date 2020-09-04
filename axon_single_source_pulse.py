from main import ExternalPotential
from parameters import monophasic_pulse_params
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
import sys

root_folder = os.path.abspath(join(os.path.dirname(__file__), '..'))
cell_models_folder = join(os.path.dirname(__file__), "cell_models")

model_path = join(cell_models_folder, 'unmyelinated_axon.hoc')

cell_parameters = {
    'morphology': model_path,
    'nsegs_method': "lambda_f",
    'lambda_f': 1000.,
    'v_init': -65,
    'passive': False,
    'dt': self.dt,  # [ms] Should be a power of 2
    'tstart': -self.cut_off,  # [ms] Simulation start time
    'tstop': self.tstop,  # [ms] Simulation end time
    "pt3d": True,
    "extracellular": True,
}

cell = LFPy.Cell(**cell_parameters)
