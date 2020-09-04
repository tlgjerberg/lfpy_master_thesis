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

    def _extra_celular_stimuli(self, elec_params):
        x0, y0, z0 = elec_params['position']
        sigma = elec_params['sigma']
        start_time = elec_params['start_time']
        stop_time = elec_params['stop_time']

        ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * sigma *
                                                      np.sqrt((x0 - x)**2 +
                                                              (y0 - y)**2 +
                                                              (z0 - z)**2)))

        n_tsteps = int(self.tstop / self.dt + 1)
        t = np.arange(n_tsteps) * self.dt
        pulse = np.zeros(n_tsteps)
        start_idx = np.argmin(np.abs(t - start_time))
        stop_idx = np.argmin(np.abs(t - stop_time))
        pulse[start_idx:stop_idx] = elec_params['amp'] * 1000  # Test * 1000

        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = ext_field(cell.xmid, cell.ymid, cell.zmid).reshape(
            cell.totnsegs, 1) * pulse.reshape(1, n_tsteps)

        cell.insert_v_ext(v_cell_ext, t)

        return ext_field, pulse

    def plot_ext_stim(self):

        plt.show()
