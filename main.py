import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
import sys
from os.path import join
from matplotlib.patches import Ellipse

"""
    Parameters: LFPy Cell object, dict containing electrode parameters

    Returns:
    """


class ExternalPotential:
    def __init__(self, cell_params):

        self.cell_params = cell_params
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cell_name = cell_params['cell_name']

    def extra_cellular_stimuli(self, cell, elec_params):

        self.cell = cell
        self.elec_params = elec_params
        self.amp = elec_params['pulse_amp']
        self.x0, self.y0, self.z0 = elec_params['positions'][0]
        sigma = elec_params['sigma']
        start_time = elec_params['start_time']
        stop_time = elec_params['stop_time']

        # Calculating external field
        self.ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * sigma *
                                                           np.sqrt((self.x0 - x)**2 +
                                                                   (self.y0 - y)**2 +
                                                                   (self.z0 - z)**2)))

        # Generating time steps of simulation
        n_tsteps = int(self.tstop / self.dt + 1)
        t = np.arange(n_tsteps) * self.dt

        # Setting pulse change at a given time interval
        self.pulse = np.zeros(n_tsteps)
        start_idx = np.argmin(np.abs(t - start_time))
        stop_idx = np.argmin(np.abs(t - stop_time))
        self.pulse[start_idx:stop_idx] = elec_params['pulse_amp'] * 1000

        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.xmid, cell.ymid, cell.zmid).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        cell.insert_v_ext(v_cell_ext, t)

        # return ext_field, pulse

    def plot_cellsim(self):

        cell_plot_idxs = [0, int(self.cell.totnsegs / 2),
                          self.cell.totnsegs - 1]
        cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
            1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}

        plt.figure(figsize=(16, 9))

        v_field_ext = np.zeros((50, 200))
        xf = np.linspace(np.min(self.cell.xend), np.max(self.cell.xend), 50)
        zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

        for xidx, x in enumerate(xf):

            for zidx, z in enumerate(zf):
                v_field_ext[xidx, zidx] = self.ext_field(x, 0, z) * self.amp

        plt.subplots_adjust(hspace=0.5)
        plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='z [$\mu m$]',
                    xlim=[-400, 400], xticks=[-400, 0, 400], title='Green dots: Measurement points')
        plt.imshow(v_field_ext.T, extent=[np.min(self.cell.xend), np.max(self.cell.xend), np.min(self.cell.zend), np.max(self.cell.zend)],
                   origin='lower', interpolation='nearest', cmap=plt.cm.bwr_r, vmin=-150, vmax=150)

        l, = plt.plot(self.x0, self.z0, 'y*', ms=12)
        plt.legend([l], ["point current source"], frameon=False)

        # Plotting the membrane potentials
        plt.subplot(222, title='Membrane potential',
                    xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
        [plt.plot(self.cell.tvec, self.cell.vmem[idx, :], c=cell_plot_colors[idx], lw=2)
         for idx in cell_plot_idxs]

        ax1 = plt.subplot(224, ylim=[-2 * np.max(np.abs(self.pulse / 1000)), 2 * np.max(np.abs(self.pulse / 1000))],
                          ylabel='$\mu$A', title='Injected current')
        ax1.plot(self.cell.tvec, self.pulse / 1000)
        plt.show()
