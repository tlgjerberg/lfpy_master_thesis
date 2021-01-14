import numpy as np
import neuron
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import LFPy
import os
import sys
import time
from os.path import join
from matplotlib.patches import Ellipse
from matplotlib.cm import ScalarMappable


"""
 transmembrane current i soma/utvalgte punkter som func av tid x
membran pot som func av tid

 snapshots ved start og stop

"""


class ExternalPotentialSim:

    def __init__(self, cell_params):

        self.root_folder = os.path.dirname(__file__)

        self.cell_params = cell_params
        self.save_folder = join(
            self.root_folder, cell_params['save_folder_name'])
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
        self.cell_name = cell_params['cell_name']
        self.cell_dist_to_top = cell_params['cell_dist_to_top']
        self.x_shift = cell_params['x_shift']
        self.z_rot = cell_params['z_rot']

    def return_cell(self, cell_models_folder, passive=True):

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
            self.v_init = cell_parameters['v_init']
            # print(self.cell.zmid[1] - self.cell.zmid[0])
            if not passive:
                neuron.h('forall insert hh')

        elif self.cell_name == 'Hallermann':
            model_path = join(cell_models_folder, 'HallermannEtAl2012')

            neuron.load_mechanisms(model_path)

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
                'custom_code': [join(model_path, 'Cell parameters.hoc'),
                                join(model_path, 'charge.hoc')]
            }
            cell = LFPy.Cell(**cell_parameters)
            self.cell = cell
            self.v_init = cell_parameters['v_init']

            # cell.set_pos(z=-np.max(cell.zend) -
            #              self.cell_dist_to_top, x=self.x_shift)
            # cell.set_rotation(z=self.z_rot)
            self.cell.set_pos(z=-self.cell_dist_to_top)
            # Default rotation
            self.cell.set_rotation(x=4.729, y=-3.05, z=-3)
            # Axon y-coordinate close to 0
            # self.cell.set_rotation(x=4.9, y=-3.166, z=-3)
            # Apical dendtrite measurement point at aprrox y=0
            # self.cell.set_rotation(x=4.788, y=-3.166, z=-3)

    def create_measure_points(self, coords):
        """
        Setting measurement of membrane potential in compartments closest to
        coordinates.
        """
        measure_pnts = []
        self.elec_positions = []

        for i in range(coords.shape[0]):
            x, y, z = coords[i, :]
            print(i, self.cell.get_closest_idx(x, y, z))
            measure_pnts.append(self.cell.get_closest_idx(x, y, z))

        self.measure_pnts = np.array(measure_pnts)

        # Automatically setting electrode at a given distance from the measurement points
        if self.cell_name == 'Hallermann':

            for mp in self.measure_pnts:
                self.elec_positions.append(
                    np.array([int(self.cell.xmid[mp]) - 50, int(self.cell.ymid[mp]),
                              int(self.cell.zmid[mp])], dtype=float))

    def _calc_point_sources_field(self, elec_params):
        pass

    def extracellular_stimuli(self, elec_params):
        """
        Parameters: LFPy Cell object, dict of electrode parameters

        Returns:
        """
        self.elec_params = elec_params
        self.amp = elec_params['pulse_amp']  # External current amplitide
        # Electrode position
        print('ex_stim', elec_params['positions'])
        self.x0, self.y0, self.z0 = elec_params['positions']
        sigma = elec_params['sigma']
        self.start_time = elec_params['start_time']
        self.stop_time = elec_params['stop_time']

        # Calculating external field function
        self.ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * sigma *
                                                           np.sqrt((self.x0 - x)**2 +
                                                                   (self.y0 - y)**2 +
                                                                   (self.z0 - z)**2)))

        # Generating time steps of simulation
        n_tsteps = int(self.tstop / self.dt + 1)
        t = np.arange(n_tsteps) * self.dt

        # Setting pulse change at a given time interval
        self.pulse = np.zeros(n_tsteps)
        self.start_idx = np.argmin(np.abs(t - self.start_time))
        self.stop_idx = np.argmin(np.abs(t - self.stop_time))
        self.pulse[self.start_idx:self.stop_idx] = elec_params['pulse_amp']

        # Applying the external field function to the cell simulation
        v_cell_ext = np.zeros((self.cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(self.cell.xmid, self.cell.ymid, self.cell.zmid).reshape(
            self.cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        self.cell.insert_v_ext(v_cell_ext, t)

    def run_cell_simulation(self):
        self.cell.simulate(rec_vmem=True, rec_imem=True)

    def _find_steady_state(self):
        """
        Detect steady state potential
        """
        self.v_ss = np.max(self.cell.vmem)
        self.v_ss_idx = np.argmax(self.cell.vmem)
        find_diff = np.diff(self.cell.vmem)

    def _dV(self):
        self._find_steady_state()
        self.dV = self.v_ss - self.v_init

    def _record_dist_to_electrode(self, coords):

        self.record_dist = np.zeros(coords.shape[0])

        for idx in range(coords.shape[0]):
            self.record_dist[idx] = np.sum(np.absolute(
                coords[idx, :] - np.array([self.x0, self.y0, self.z0])))

        return self.record_dist

    def find_time_constant(self):
        pass

    def run_ext_sim(self, cell_models_folder, elec_params, current_amps, positions, coords, stop_time, passive=False):

        self.return_cell(cell_models_folder, passive)
        self.create_measure_points(coords)
        elec_dists = np.zeros((len(self.elec_positions), coords.shape[0]))
        ss_pot = np.zeros(len(self.elec_positions))
        dV = np.zeros(len(self.elec_positions))
        # Neuron activation after cell object has been created

        for I in current_amps:

            elec_params['pulse_amp'] = I

            for idx, pos in enumerate(self.elec_positions):
                print('idx', idx)
                elec_params['positions'] = pos
                self.extracellular_stimuli(elec_params)
                self.run_cell_simulation()
                self.plot_cellsim()
                self._find_steady_state()
                ss_pot[idx] = self.v_ss
                self._dV()
                dV[idx] = self.dV
                elec_dists[idx] = self._record_dist_to_electrode(coords)

        self.plot_steady_state(elec_dists[:, 0], ss_pot)
        self.plot_dV(elec_dists[:, 0], dV)
        self.create_measure_points(coords)

        # Freeing up some variables
        I = None
        pos = None

    def run_current_sim(self, cell_models_folder, elec_params, current_amps, positions, stop_time, passive=False):
        self.return_cell(cell_models_folder)

        # Neuron activation after cell object has been created
        if not passive:
            neuron.h('forall insert hh')

        elec_params['positions'] = positions
        elec_params['pulse_amp'] = current_amps

        self.extracellular_stimuli(elec_params)
        self.run_cell_simulation()
        self.plot_currents()

    def plot_morphology(self):

        # Defining figure frame and parameters
        self.fig.subplots_adjust(hspace=0.5, left=0.5, wspace=0.5, right=0.96,
                                 top=0.9, bottom=0.1)

        # Adding axes with appropriate parameters
        self.ax_m = self.fig.add_axes(self.morph_ax_params, aspect=1, frameon=False,
                                      xticks=[], yticks=[], ylim=[-700, 1100], xlim=[-300, 300])

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

            self.ax_m.plot([self.cell.xstart[idx], self.cell.xend[idx]],
                           [self.cell.zstart[idx], self.cell.zend[idx]], '-',
                           c=c, clip_on=True, lw=np.sqrt(self.cell.diam[idx]) * 1)

        lines = []
        for name in used_clrs:
            l, = self.ax_m.plot([0], [0], lw=2, c=sec_clrs[name])
            lines.append(l)
        self.ax_m.legend(lines, used_clrs, frameon=False,
                         fontsize=8, loc=(0.05, 0.0), ncol=2)

        # Plotting dots at the middle of a given section in its given color
        [self.ax_m.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o',
                        c=self.cell_plot_colors[idx], ms=13) for idx in self.cell_plot_idxs]

        self.ax_m.text(20, 40, "Cortical electrode\n(R={} $\mu$m)".format(self.elec_params["electrode_radii"]),
                       fontsize=9, ha='center')

        # for e_idx in range(len(self.elec_params["positions"])):
        #     ellipse_pos = [self.elec_params["positions"][e_idx]
        #                    [0], self.elec_params["positions"][e_idx][2]]
        # print(self.elec_params["positions"])
        ellipse_pos = [self.elec_params["positions"]
                       [0], self.elec_params["positions"][2]]

        self.ax_m.add_artist(Ellipse(ellipse_pos, width=2 * self.elec_params["electrode_radii"],
                                     height=self.elec_params["electrode_radii"] / 5, fc='gray', ec='black'))

    def plot_external_field(self):
        # Adding external field visualization to cell morphology figure
        v_field_ext = np.zeros((200, 200))
        x = np.linspace(-500, 500, 200)
        z = np.linspace(-500, 1000, 200)
        # x = np.linspace(np.min(self.cell.xend), np.max(self.cell.xend), 50)
        # z = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)
        xf, zf = np.meshgrid(x, z)

        for xidx, xi in enumerate(x):

            for zidx, zi in enumerate(z):

                v_field_ext[xidx, zidx] = self.ext_field(xi, 0, zi) * self.amp

        vmax = np.max(np.abs(v_field_ext)) / 5
        # ax_cb = plt.gca()
        # im_p = ax_cb.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
        #                     origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)
        #
        divider = make_axes_locatable(self.ax_m)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        # plt.colorbar(im_p, cax=cax, label='mV')
        ext_field_im = self.ax_m.pcolormesh(
            xf, zf, v_field_ext.T, cmap='bwr', vmin=-vmax, vmax=vmax, shading='auto')
        self.fig.colorbar(ScalarMappable(
            norm=None, cmap='bwr'), cax=cax, ax=self.ax_m)
        # plt.show()

    def plot_membrane_potential(self, placement):

        ax_vm = self.fig.add_axes(placement,  # ylim=[-120, 50],
                                  xlim=[0, self.tstop], xlabel="Time (ms)")

        ax_vm.set_ylabel("Membrane\npotential (mV)", labelpad=-3)

        # if type(self.spike_time_idxs) == int:
        #     ax_vm.axvline(self.cell.tvec[self.spike_time_idxs], c='r', ls='--')

        # mark_subplots([ax_stim, ax_vm], "BC", xpos=-0.02, ypos=0.98)
        [ax_vm.plot(self.cell.tvec, self.cell.vmem[idx],
                    c=self.cell_plot_colors[idx], lw=0.5) for idx in self.cell_plot_idxs]

    def plot_current_pulse(self, placement):

        ax_stim = self.fig.add_axes(placement, xlim=[0, self.tstop],
                                    ylabel="Stimuli\ncurrent ($\mu$A)", xlabel="Time (ms)")
        # ax_stim.set_ylabel("$\mu$A", labelpad=-2)
        ax_stim.plot(self.cell.tvec, self.pulse / 1000, lw=0.5)

    def plot_cellsim(self):
        # Simulating cell after all parameters and field has been added
        self.fig = plt.figure(figsize=[10, 8])

        # for m in self.measure_pnts:
        #     print((self.cell.xmid[m], self.cell.ymid[m], self.cell.zmid[m]))

        self.cell_plot_idxs = self.measure_pnts.astype(
            dtype='int')  # List of measurement points

        self.cell_plot_colors = idx_clr = {idx: [
            'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        # Setting size and location of plotted morphology
        self.morph_ax_params = [0.1, 0.05, 0.2, 0.90]

        self.plot_morphology()
        self.plot_external_field()

        # # Setting size and location of plotted potentials and current
        ax_top = 0.90
        ax_h = 0.30
        ax_w = 0.6
        ax_left = 0.3
        stim_axes_placement = [ax_left, ax_top - ax_h, ax_w, ax_h]
        mem_axes_placement = [ax_left, ax_top - ax_h - 0.47, ax_w, ax_h]

        self.plot_membrane_potential(mem_axes_placement)
        self.plot_current_pulse(stim_axes_placement)

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        self.fig.savefig(join(
            self.save_folder, f'ext_field_point_amp={self.amp}uA_x={self.x0}_z={self.z0}.png'))
        # plt.show()
        plt.close(fig=self.fig)

    def plot_steady_state(self, elec_abs_dists, steady_state):

        fig, ax = plt.subplots()
        ax.set_xlabel('Electrode Distance from Cell Origin ($\mu m$)')
        ax.set_ylabel('Steady State Potential (mV)')
        ax.plot(elec_abs_dists, steady_state, '-o')

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        plt.savefig(
            join(self.save_folder, 'potential_electrode_distance.png'), dpi=300)
        # plt.show()

    def plot_dV(self, elec_dists, dV):

        # elec_dists = np.log(elec_dists)
        # dV = np.log(dV)

        fig, ax = plt.subplots()
        ax.set_xlabel('Electrode Distance from Cell Origin ($\mu m$)')
        ax.set_ylabel('dV (mV)')
        ax.loglog(elec_dists, dV, '-o')

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        fig.savefig(
            join(self.save_folder, 'dV_electrode_distance.png'), dpi=300)

    def plot_currents(self):

        # Midpoint index for pulse as an extra test point in time
        self.mid_idx = (self.stop_idx + self.start_idx) // 2

        # Time indices for snapshots currents
        timepoints = np.array(
            [self.start_idx + 1, self.mid_idx, self.stop_idx])

        # Extracting axial currents and their coordinates
        ax_current, _, pos_coord = self.cell.get_axial_currents_from_vmem(
            timepoints=timepoints)

        save_folder = 'axon_bisc_currents'

        self.save_folder = join(self.root_folder, save_folder)

        self.fig = plt.figure()  # figsize=[18, 8]

        self.cell_plot_idxs = self.measure_pnts.astype(
            dtype='int')  # List of measurement points

        self.cell_plot_colors = idx_clr = {idx: [
            'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        self.morph_ax_params = [0.6, 0.05, 0.2, 0.90]

        self.plot_morphology()
        self.plot_external_field()

        ax_top = 0.90
        ax_h = 0.30
        ax_w = 0.3
        ax_left = 0.2
        mem_axes_placement = [ax_left, ax_top - ax_h - 0.47, ax_w, ax_h]

        ax_tmc = self.fig.add_axes([ax_left, ax_top - ax_h, ax_w, ax_h],  # ylim=[-120, 50],
                                   xlim=[0, self.tstop], ylabel='Transmembrane Current (mV)', xlabel="Time (ms)")

        # [ax_tmc.plot(self.cell.tvec, self.cell.imem[idx, :],
        #              c=self.cell_plot_colors[idx], lw=0.5) for idx in self.cell_plot_idxs]

        self.plot_membrane_potential(mem_axes_placement)

        # Plotting snapshots at start and stop times of current pulse
        fig_snap1, ax_snap1 = plt.subplots()
        ax_snap1.plot(
            self.cell.imem[:, self.start_idx + 1], self.cell.zmid, 'o-')
        ax_snap1.axvline(0, ls="--", c='grey')
        ax_snap1.set_xlabel(f'Current at time {self.start_idx + 1} (nA)')
        ax_snap1.set_ylabel('Cell Compartments in z direction')

        fig_snap2, ax_snap2 = plt.subplots()
        ax_snap2.axvline(0, ls="--", c='grey')
        ax_snap2.plot(self.cell.imem[:, self.stop_idx], self.cell.zmid, 'o-')
        ax_snap1.set_xlabel(f'Current at time {self.stop_idx} (nA)')
        ax_snap1.set_ylabel('Cell Compartments in z direction')

        # Plotting axial current along the z-direction of
        fig_axial, ax_axial = plt.subplots()
        ax_axial.plot(ax_current[:, 0], pos_coord[:, 2])
        ax_axial.axvline(0, ls="--", c='grey')
        ax_axial.set_xlim([-6, 6])
        ax_axial.set_xlabel(f'Axial Current at time {timepoints[0]} ms (nA)')
        ax_axial.set_ylabel('Cell Compartments Along Axon (z direction)')

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        self.fig.savefig(join(
            self.save_folder, f'transmembrane_current_amp={self.amp}.svg'), dpi=300)
        fig_snap1.savefig(join(
            self.save_folder, f'transmembrane_current_snapshot t={self.start_idx + 1}.png'), dpi=300)
        fig_snap2.savefig(join(
            self.save_folder, f'transmembrane_current_snapshot t={self.stop_idx}.png'), dpi=300)
        fig_axial.savefig(join(
            self.save_folder, f'axial_current_soma_snapshot t={timepoints[0]}.png'), dpi=300)

        # plt.show()

        plt.close(fig=self.fig)
