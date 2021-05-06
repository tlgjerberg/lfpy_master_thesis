from os.path import join
import os
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from main import ExternalPotentialSim
import matplotlib
from mpl_toolkits import mplot3d
# matplotlib.use("AGG")


class PlotSimulations(ExternalPotentialSim):
    def __init__(self, cell_params, elec_params, cell_vmem=None, cell_tvec=None):
        super().__init__(cell_params, elec_params)
        # self.morph_ax_params = morph_ax_params  # Move to plot_morphology?
        self.vmem = cell_vmem
        self.tvec = cell_tvec

    def plot_morphology(self, cell, morph_ax_params):

        # Adding axes with appropriate parameters
        self.ax_m = self.fig.add_axes(morph_ax_params, aspect=1, frameon=False,
                                      xticks=[], yticks=[], ylim=self.ylim, xlim=self.xlim)
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
        for idx in range(cell.totnsegs):
            sec_name = cell.get_idx_name(idx)[1]
            # print(sec_name)
            # c = 'k'
            for ax_name in possible_names:
                if ax_name in sec_name:
                    # print(ax_name, sec_name)
                    c = sec_clrs[ax_name]
                    if not ax_name in used_clrs:
                        used_clrs.append(ax_name)

            # Plotting the neuron morphology with colored lines matching segments
            self.ax_m.plot([cell.x[idx][0], cell.x[idx][1]],
                           [cell.z[idx][0], cell.z[idx][1]], '-',
                           c=c, clip_on=True, lw=2)

        # Adding discriptors below the morphology for each kind of segement
        lines = []
        for name in used_clrs:
            l, = self.ax_m.plot([0], [0], lw=2, c=sec_clrs[name])
            lines.append(l)
        self.ax_m.legend(lines, used_clrs, frameon=False,
                         fontsize=8, loc=(0.05, 0.0), ncol=2)

        # Plotting dots at the middle of a given section in its given color
        [self.ax_m.plot(cell.x.mean(axis=1)[idx], cell.z.mean(axis=1)[idx], 'o',
                        c=self.cell_plot_colors[idx], ms=7) for idx in self.cell_plot_idxs]

    def morphology_3D(self, cell):

        # List of measurement points
        self.cell_plot_idxs = self.measure_pnts.astype(dtype='int')

        self.cell_plot_colors = {idx: [
            'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        fig_3dmorph = plt.figure()
        ax_3dmorph = plt.axes(projection="3d")

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

        for idx in range(cell.totnsegs):
            sec_name = cell.get_idx_name(idx)[1]
            # print(sec_name)
            # c = 'k'
            for ax_name in possible_names:
                if ax_name in sec_name:
                    # print(ax_name, sec_name)
                    c = sec_clrs[ax_name]
                    if not ax_name in used_clrs:
                        used_clrs.append(ax_name)

            ax_3dmorph.plot([cell.x[idx][0], cell.x[idx][1]],
                            [cell.y[idx][0], cell.y[idx][1]],
                            [cell.z[idx][0], cell.z[idx][1]], '-',
                            c=c, clip_on=True, lw=2)
        lines = []
        for name in used_clrs:
            l, = ax_3dmorph.plot([0], [0], lw=2, c=sec_clrs[name])
            lines.append(l)
        ax_3dmorph.legend(lines, used_clrs, frameon=False,
                          fontsize=8, loc=(0.05, 0.0), ncol=2)

        ax_3dmorph.legend(lines, used_clrs, frameon=False,
                          fontsize=8, loc=(0.05, 0.0), ncol=2)
        ax_3dmorph.set_xlabel(r'x[$\mu m$]')
        ax_3dmorph.set_ylabel(r'y[$\mu m$]')
        ax_3dmorph.set_zlabel(r'z[$\mu m$]')

        # ellipse_pos = [self.elec_params["positions"][0],
        #                self.elec_params["positions"][1],
        #                self.elec_params["positions"][2]]
        # # ax_3dmorph.plot(ellipse_pos[0], ellipse_pos[1], ellipse_pos[2], 'o')
        [ax_3dmorph.plot(cell.x.mean(axis=1)[idx], cell.y.mean(axis=1)[idx], cell.z.mean(axis=1)[idx], 'o',
                         c=self.cell_plot_colors[idx], ms=3) for idx in self.cell_plot_idxs]

        fig_3dmorph.savefig(
            join(self.save_folder, f'{self.cell_name}_{self.amp}mA_elec_pos={self.elec_pos}.png'))

    def draw_electrode(self):
        ellipse_pos = [self.elec_params["positions"]
                       [0], self.elec_params["positions"][2]]

        self.ax_m.add_artist(Ellipse(ellipse_pos, width=10 * self.elec_params["electrode_radii"],
                                     height=self.elec_params["electrode_radii"] / 5, fc='gray', ec='black'))

    def plot_external_field(self, cell):

        # Adding external field visualization to cell morphology figure
        field_x_dim = abs(self.xlim[0]) + abs(self.xlim[1])
        field_y_dim = abs(self.ylim[0]) + abs(self.ylim[1])
        v_field_ext = np.zeros((field_x_dim, field_y_dim))
        x = np.linspace(self.xlim[0], self.xlim[1], field_x_dim)
        z = np.linspace(self.ylim[0], self.ylim[1], field_y_dim)
        # x = np.linspace(np.min(cell.xend), np.max(cell.xend), 50)
        # z = np.linspace(np.min(cell.zend), np.max(cell.zend), 200)
        xf, zf = np.meshgrid(x, z)

        for xidx, xi in enumerate(x):

            for zidx, zi in enumerate(z):

                v_field_ext[xidx, zidx] = self.ext_field(xi, 0, zi) * self.amp

        vmax = np.max(np.abs(v_field_ext)) / 5

        # Creating a field image
        self.ext_field_im = self.ax_m.pcolormesh(
            xf, zf, v_field_ext.T, cmap='bwr', vmin=-vmax, vmax=vmax, shading='auto')

        self.add_colorbar()

    def add_colorbar(self):

        # Divide axes of field and colorbar
        divider = make_axes_locatable(self.ax_m)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        # Add colorbar corresonding to field strength in mV
        self.fig.colorbar(self.ext_field_im, pad=0.05, ax=self.ax_m, cax=cax)

    def plot_membrane_potential(self, placement):

        ax_vm = self.fig.add_axes(placement,  # ylim=[-120, 50],
                                  xlim=[0, self.tstop], xlabel="Time (ms)")

        ax_vm.set_ylabel("Membrane\npotential (mV)", labelpad=-3)

        # mark_subplots([ax_stim, ax_vm], "BC", xpos=-0.02, ypos=0.98)
        [ax_vm.plot(self.tvec, self.vmem[idx],
                    c=self.cell_plot_colors[idx], lw=0.5) for idx in self.cell_plot_idxs]

    def plot_current_pulse(self, placement):

        ax_stim = self.fig.add_axes(placement, xlim=[0, self.tstop],
                                    ylabel="Stimuli\ncurrent ($\mu$A)", xlabel="Time (ms)")
        # ax_stim.set_ylabel("$\mu$A", labelpad=-2)
        ax_stim.plot(self.tvec, self.pulse / 1000, lw=0.5)

    def plot_cellsim(self, cell, com_coords, z_rot, morph_ax_params, xlim=[-500, 500], ylim=[-300, 1200]):
        """
        To add:
        if placement:
            plot in one figure
        else:
            separate figures?
        """
        self.xlim = xlim
        self.ylim = ylim

        # Recreating the cell object used in simulatios without running simul
        cell = self.return_cell(self.cell_models_folder)
        self.create_measure_points(cell, com_coords)
        self.extracellular_stimuli(cell)

        # Simulating cell after all parameters and field has been added
        self.fig = plt.figure(figsize=[10, 8])

        # Defining figure frame and parameters for combined figure
        self.fig.subplots_adjust(hspace=0.5, left=0.5, wspace=0.5, right=0.96,
                                 top=0.9, bottom=0.1)

        self.cell_plot_idxs = self.measure_pnts.astype(
            dtype='int')  # List of measurement points

        self.cell_plot_colors = {idx: [
            'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        self.plot_morphology(cell, morph_ax_params)
        self.plot_external_field(cell)

        # # Setting size and location of plotted potentials and current
        ax_top = 0.90
        ax_h = 0.30
        ax_w = 0.45
        ax_left = 0.5
        stim_axes_placement = [ax_left, ax_top - ax_h, ax_w, ax_h]
        mem_axes_placement = [ax_left, ax_top - ax_h - 0.47, ax_w, ax_h]

        self.plot_membrane_potential(mem_axes_placement)
        self.plot_current_pulse(stim_axes_placement)
        self.draw_electrode()

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        self.fig.savefig(join(
            self.save_folder, f'point_source_{self.cell_name}_z_rot={self.z_rot}_point_amp={self.amp}uA_x={self.x0}_z={self.z0}.png'))
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

    def plot_currents(self, cell):

        # Midpoint index for pulse as an extra test point in time
        self.mid_idx = (self.stop_idx + self.start_idx) // 2

        # Time indices for snapshots currents
        timepoints = np.array(
            [self.start_idx + 1, self.mid_idx, self.stop_idx])

        # Extracting axial currents and their coordinates
        ax_current, _, pos_coord = cell.get_axial_currents_from_vmem(
            timepoints=timepoints)

        save_folder = 'axon_bisc_currents'

        self.save_folder = join(self.root_folder, save_folder)

        self.fig = plt.figure()  # figsize=[18, 8]

        # List of measurement points
        self.cell_plot_idxs = self.measure_pnts.astype(dtype='int')

        self.cell_plot_colors = {idx: [
            'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        # [0.05, 0.05, 0.2, 0.90]
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

        # [ax_tmc.plot(cell.tvec, cell.imem[idx, :],
        #              c=cell_plot_colors[idx], lw=0.5) for idx in self.cell_plot_idxs]

        self.plot_membrane_potential(mem_axes_placement)

        # Plotting snapshots at start and stop times of current pulse
        fig_snap1, ax_snap1 = plt.subplots()
        ax_snap1.plot(
            cell.imem[:, self.start_idx + 1], cell.zmid, 'o-')
        ax_snap1.axvline(0, ls="--", c='grey')
        ax_snap1.set_xlabel(f'Current at time {self.start_idx + 1} (nA)')
        ax_snap1.set_ylabel('Cell Compartments in z direction')

        fig_snap2, ax_snap2 = plt.subplots()
        ax_snap2.axvline(0, ls="--", c='grey')
        ax_snap2.plot(cell.imem[:, self.stop_idx], cell.zmid, 'o-')
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
            self.save_folder, f'transmembrane_current_amp={self.amp}.png'), dpi=300)
        fig_snap1.savefig(join(
            self.save_folder, f'transmembrane_current_snapshot t={self.start_idx + 1}.png'), dpi=300)
        fig_snap2.savefig(join(
            self.save_folder, f'transmembrane_current_snapshot t={self.stop_idx}.png'), dpi=300)
        fig_axial.savefig(join(
            self.save_folder, f'axial_current_soma_snapshot t={timepoints[0]}.png'), dpi=300)

        # plt.show()

        plt.close(fig=self.fig)

    def plot_double_morphology(self, cell_models_folder, z_rot, msre_coords):

        morph_ax_params1 = [0.05, 0.1, 0.5, 0.9]
        morph_ax_params2 = [0.48, 0.1, 0.5, 0.9]

        self.fig = plt.figure()
        self.z_rot = z_rot[0]
        cell1 = self.return_cell(cell_models_folder)
        self.create_measure_points(cell1, msre_coords)

        self.cell_plot_idxs = self.measure_pnts.astype(dtype='int')
        self.cell_plot_colors = {idx: [
            'green', 'cyan', 'orange', 'b', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}

        self.plot_morphology(cell1, morph_ax_params1)
        self.extracellular_stimuli(cell1)

        if self.elec_pos[0] > 0:
            self.draw_electrode()
        cell1.__del__()

        self.z_rot = z_rot[1]
        cell2 = self.return_cell(cell_models_folder)
        self.cell_plot_idxs = self.measure_pnts.astype(dtype='int')
        self.cell_plot_colors = {idx: [
            'orange', 'cyan', 'b', 'green', 'purple'][num] for num, idx in enumerate(self.cell_plot_idxs)}
        self.create_measure_points(cell2, msre_coords)
        self.plot_morphology(cell2, morph_ax_params2)

        if self.elec_pos[0] < 0:
            self.draw_electrode()
        cell2.__del__()

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)
        # plt.show()
        self.fig.savefig(join(
            self.save_folder, f'point_source_{self.cell_name}_z_rot={self.z_rot}_point_amp={self.amp}uA_elec_pos={self.elec_pos}.png'))

        plt.close(fig=self.fig)

    def plot_double_mem_pot(self, cell_vmem, cell_tvec):
        mem_axes_placement1 = [0.2, 0.6, 0.6, 0.3]
        mem_axes_placement2 = [0.2, 0.2, 0.6, 0.3]
        self.fig = plt.figure()

        self.plot_membrane_potential(mem_axes_placement1)

        self.tvec = cell_tvec
        self.vmem = cell_vmem

        self.plot_membrane_potential(mem_axes_placement2)

        self.fig.savefig(join(
            self.save_folder, f'double_vmem_{self.cell_name}_z_rot={self.z_rot}_point_amp={self.amp}uA_elec_pos={self.elec_pos}.png'))
        plt.show()
        plt.close(fig=self.fig)
