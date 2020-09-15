import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
import sys
from os.path import join
from matplotlib.patches import Ellipse


class ExternalPotentialSim:

    def __init__(self, cell_params):

        self.cell_params = cell_params
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
        self.cell_name = cell_params['cell_name']
        self.cell_dist_to_top = cell_params['cell_dist_to_top']
        self.x_shift = cell_params['x_shift']
        self.z_rot = cell_params['z_rot']

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
            cell = LFPy.Cell(**cell_parameters)
            self.cell = cell

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
            cell = LFPy.Cell(**cell_parameters)
            self.cell = cell

            # cell.set_pos(z=-np.max(cell.zend) -
            #              self.cell_dist_to_top, x=self.x_shift)
            # cell.set_rotation(z=self.z_rot)
            self.cell.set_pos(z=-self.cell_dist_to_top)
            self.cell.set_rotation(x=4.729, y=-3.166, z=-3)

        # return cell

    def extra_cellular_stimuli(self, elec_params):
        """
        Parameters: LFPy Cell object, dict of electrode parameters

        Returns:
        """
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

        v_cell_ext = np.zeros((self.cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(self.cell.xmid, self.cell.ymid, self.cell.zmid).reshape(
            self.cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        self.cell.insert_v_ext(v_cell_ext, t)

        # return ext_field, pulse

    def plot_cellsim(self):
        self.cell.simulate(rec_vmem=True)

        # Setting cell compartments to measure AP
        cell_plot_idxs = [0, int(self.cell.totnsegs / 2),
                          self.cell.totnsegs - 1]
        cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
            1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}

        plt.figure(figsize=(16, 9))

        v_field_ext = np.zeros((50, 200))
        xf = np.linspace(-500, 500, 50)
        zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

        for xidx, x in enumerate(xf):

            for zidx, z in enumerate(zf):
                v_field_ext[xidx, zidx] = self.ext_field(x, 0, z) * self.amp

        vmax = np.max(np.abs(v_field_ext)) / 5
        plt.subplots_adjust(hspace=0.5)
        plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]',
                    xlim=[-500, 500], xticks=[-500, 0, 500], title='Green dots: Measurement points')
        plt.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
                   origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)

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

    def plot_cellsim_alt(self, measure_idxs):
        self.cell.simulate(rec_vmem=True)

        cell_plot_idxs = measure_idxs.astype(
            dtype='int')  # List of measurement points
        cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
            1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}

        # Defining figure frame and parameters
        fig = plt.figure(figsize=[18, 8])
        fig.subplots_adjust(hspace=0.5, left=0.0, wspace=0.4, right=0.96,
                            top=0.97, bottom=0.1)

        # Adding axes with appropriate parameters
        ax_m = fig.add_axes([-0.01, 0.05, 0.2, 0.97], aspect=1, frameon=False,
                            xticks=[], yticks=[], ylim=[-1900, 300], xlim=[-300, 300])

        # Names of different neuron parts and color codings for each
        possible_names = ["Myelin", "axon", "Unmyelin", "Node", "hilloc",
                          "hill", "apic", "dend", "soma"]
        sec_clrs = {"Myelin": 'olive',
                    "dend": '0.3',
                    "soma": 'k',
                    'apic': '0.6',
                    "axon": 'lightgreen',
                    "Unmyelin": 'salmon',
                    "Node": 'r',
                    "hilloc": 'lightblue',
                    "hill": 'pink', }
        used_clrs = []

        # PLOTTING CELL MORPHOLOGY

        # Sets each segment to the color matching the name set by sec_clrs
        for idx in range(self.cell.totnsegs):
            sec_name = self.cell.get_idx_name(idx)[1]
            # print(sec_name)
            # c = 'k'
            for ax_name in possible_names:
                if ax_name in sec_name:
                    # print(ax_name, sec_name)
                    c = sec_clrs[ax_name]
                    if not ax_name in used_clrs:
                        used_clrs.append(ax_name)

            ax_m.plot([self.cell.xstart[idx], self.cell.xend[idx]],
                      [self.cell.zstart[idx], self.cell.zend[idx]], '-',
                      c=c, clip_on=True, lw=np.sqrt(self.cell.diam[idx]) * 1)

        lines = []
        for name in used_clrs:
            l, = ax_m.plot([0], [0], lw=2, c=sec_clrs[name])
            lines.append(l)
        ax_m.legend(lines, used_clrs, frameon=False,
                    fontsize=8, loc=(0.05, 0.0), ncol=2)

        # Plotting dots at the middle of a given section in its given color
        [ax_m.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o',
                   c=cell_plot_colors[idx], ms=13) for idx in cell_plot_idxs]

        ax_m.text(20, 40, "Cortical electrode\n(R={} $\mu$m)".format(self.elec_params["electrode_radii"]),
                  fontsize=9, ha='center')

        for e_idx in range(len(self.elec_params["positions"])):
            ellipse_pos = [self.elec_params["positions"][e_idx]
                           [0], self.elec_params["positions"][e_idx][2]]

            ax_m.add_artist(Ellipse(ellipse_pos, width=2 * self.elec_params["electrode_radii"],
                                    height=self.elec_params["electrode_radii"] / 5, fc='gray', ec='black'))

        # Adding external field visualization to cell morphology figure
        v_field_ext = np.zeros((50, 200))
        xf = np.linspace(-500, 500, 50)
        zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

        for xidx, x in enumerate(xf):

            for zidx, z in enumerate(zf):
                v_field_ext[xidx, zidx] = self.ext_field(x, 0, z) * self.amp

        vmax = np.max(np.abs(v_field_ext)) / 5
        plt.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
                   origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)

        plt.colorbar(label='mV')
        # [plt.plot([self.cell.xstart[idx], self.cell.xend[idx]], [self.cell.zstart[idx], self.cell.zend[idx]], c='gray', zorder=1)
        #  for idx in range(self.cell.totnsegs)]
        # [plt.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o', c=cell_plot_colors[idx], ms=12)
        #  for idx in cell_plot_idxs]

        ax_top = 0.95
        ax_h = 0.3
        ax_w = 0.6
        ax_left = 0.3

        ax_vm = fig.add_axes([ax_left, ax_top - ax_h - 0.47, ax_w, ax_h],  # ylim=[-120, 50],
                             xlim=[0, self.tstop], xlabel="Time (ms)")

        ax_vm.set_ylabel("Membrane\npotential (mV)", labelpad=-3)

        # if type(self.spike_time_idxs) == int:
        #     ax_vm.axvline(self.cell.tvec[self.spike_time_idxs], c='r', ls='--')
        ax_stim = fig.add_axes([ax_left, ax_top - ax_h, ax_w, ax_h], xlim=[0, self.tstop],
                               ylabel="Stimuli\ncurrent ($\mu$A)", xlabel="Time (ms)")
        # ax_stim.set_ylabel("$\mu$A", labelpad=-2)
        ax_stim.plot(self.cell.tvec, self.pulse / 1000, lw=0.5)

        # mark_subplots([ax_stim, ax_vm], "BC", xpos=-0.02, ypos=0.98)
        [ax_vm.plot(self.cell.tvec, self.cell.vmem[idx],
                    c=cell_plot_colors[idx], lw=0.5) for idx in cell_plot_idxs]

        plt.show()
