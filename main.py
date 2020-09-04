import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
import sys


class ExternalPotential:
    def __init__(self, cell_params):

        self.cell_params = cell_params
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cell_name = cell_params['cell_name']
        self.

    def _extra_celular_stimuli(self, elec_params):
        x0, y0, z0 = elec_params['position']
        sigma = elec_params['sigma']

        ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * sigma *
                                                      np.sqrt((x0 - x)**2 +
                                                              (y0 - y)**2 +
                                                              (z0 - z)**2)))
