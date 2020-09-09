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


class ExternalPotentialSim:
    def __init__(self, cell_params):

        self.cell_params = cell_params
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
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
        self.pulse[start_idx:stop_idx] = elec_params['pulse_amp']

        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.xmid, cell.ymid, cell.zmid).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        cell.insert_v_ext(v_cell_ext, t)

        # return ext_field, pulse

    def return_cell(self, cell_models_folder):

        if self.cell_name == 'axon':
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
        elif self.cell_name == 'Hallermann':
            model_path = join(cell_models_folder, 'HallermannEtAl2012')

            neuron.load_mechanisms(cell_models_folder)

            cell_parameters = {          # various cell parameters,
                # 'morphology' : 'patdemo/cells/j4a.hoc', # Mainen&Sejnowski, 1996
                # Mainen&Sejnowski, 1996
                'morphology': join(model_path, '28_04_10_num19.hoc'),
                # 'morphology' : join('morphologies', 'axon.hoc'), # Mainen&Sejnowski, 1996
                # 'rm' : 30000.,      # membrane resistance
                # 'cm' : 1.0,         # membrane capacitance
                # 'Ra' : 150,        # axial resistance
                # 'passive_parameters':dict(g_pas=1/30., e_pas=-65),
                'v_init': -80.,    # initial crossmembrane potential
                # 'e_pas' : -65.,     # reversal potential passive mechs
                'passive': False,   # switch on passive mechs
                'nsegs_method': 'lambda_f',
                'lambda_f': 100.,
                'dt': self.dt,   # [ms] dt's should be in powers of 2 for both,
                'tstart': -self.cut_off,    # start time of simulation, recorders start at t=0
                'tstop': self.tstop,   # stop simulation at 200 ms. These can be overridden
                                    # by setting these arguments i cell.simulation()
                "extracellular": True,
                "pt3d": True,
                # 'custom_code': [join(cell_models_folder, 'Cell parameters.hoc'),
                #                 join(cell_models_folder, 'charge_only_unmyelinated.hoc')]
            }

        return cell_parameters

    def plot_cellsim(self):

        cell_plot_idxs = [0, int(self.cell.totnsegs / 2),
                          self.cell.totnsegs - 1]
        cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
            1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}

        plt.figure(figsize=(16, 9))

        v_field_ext = np.zeros((50, 200))
        xf = np.linspace(-200, 200, 50)
        zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

        for xidx, x in enumerate(xf):

            for zidx, z in enumerate(zf):
                v_field_ext[xidx, zidx] = self.ext_field(x, 0, z) * self.amp

        vmax = np.max(np.abs(v_field_ext)) / 5
        print(vmax)
        plt.subplots_adjust(hspace=0.5)
        plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]',
                    xlim=[-500, 500], xticks=[-500, 0, 500], title='Green dots: Measurement points')
        plt.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
                   origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)

        # plt.subplot(121)
        # plt.imshow(np.abs(v_field_ext.T), extent=[np.min(xf), np.max(
        #     xf), np.min(zf), np.max(zf)])

        plt.colorbar(label='mV')
        [plt.plot([self.cell.xstart[idx], self.cell.xend[idx]], [self.cell.zstart[idx], self.cell.zend[idx]], c='gray', zorder=1)
         for idx in range(self.cell.totnsegs)]
        [plt.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o', c=cell_plot_colors[idx], ms=12)
         for idx in cell_plot_idxs]

        l, = plt.plot(self.x0, self.z0, 'y*', ms=2)
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
