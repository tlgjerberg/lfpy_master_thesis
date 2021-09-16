import numpy as np
import neuron
import matplotlib
import LFPy
import os
import sys
import time
from os.path import join
from mpi4py import MPI
from neurosim import NeuronSimulation
import matplotlib.pyplot as plt


class ExternalPotentialSimulation(NeuronSimulation):

    def __init__(self, cell_params, elec_params):
        super().__init__(cell_params)

        self.amp = elec_params['pulse_amp']
        self.start_time = elec_params['start_time']
        self.stop_time = elec_params['stop_time']
        self.sigma = elec_params['sigma']  # Extracellular conductivity
        # Electrode position
        self.x0, self.y0, self.z0 = elec_params['positions']
        self.pulse_type = elec_params['pulse_type']
        self.electrode_radii = elec_params['electrode_radii']

        self._sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_shift={self.cell_dist_to_top}_z_rot={self.z_rot:.2f}_y_rot={self.y_rot:.2f}_elec_pos={self.x0}_{self.y0}_{self.z0}_t={self.stop_time}_Amp={self.amp}mA'

    def return_sim_name(self):
        """Returns a string containing simulation paramters"""

        return self._sim_name

    def update_sim_name(self):

        self._sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_shift={self.cell_dist_to_top}_z_rot={self.z_rot:.2f}_y_rot={self.y_rot:.2f}_elec_pos={self.x0}_{self.y0}_{self.z0}_t={self.stop_time}_Amp={self.amp}mA'

    def _monophaic_pulse(self, n_tsteps, t):
        # Setting monophasic pulse change at a given time interval
        self.pulse = np.zeros(n_tsteps)
        self.start_idx = np.argmin(np.abs(t - self.start_time))
        self.stop_idx = np.argmin(np.abs(t - self.stop_time))
        self.pulse[self.start_idx:self.stop_idx] = self.amp

    def insert_pulse(self, n_tsteps, t):
        if self.pulse_type == 'monophasic':
            self._monophaic_pulse(n_tsteps, t)

    def extracellular_stimuli(self, cell):
        """
        Computes the field produced by an extracellular electrode pulse and
        applies it to the cell simulation.

        Parameters: LFPy Cell object, dict of electrode parameters

        Adds an external potential to the cell object defined by a pulse
        """

        # Creating external field function
        self.ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * self.sigma *
                                                           np.sqrt((self.x0 - x)**2 +
                                                                   (self.y0 - y)**2 +
                                                                   (self.z0 - z)**2)))

        # Generating time steps for pulse in simulation
        n_tsteps = int(self.tstop / self.dt + 1)
        t = np.arange(n_tsteps) * self.dt

        self.insert_pulse(n_tsteps, t)

        # Applying the external field function to the cell simulation
        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.x.mean(axis=1), cell.y.mean(axis=1),
                                          cell.z.mean(axis=1)).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        # Adding the external potential to cell simulation
        cell.insert_v_ext(v_cell_ext, t)

    def find_steady_state_pot(self, cell_vmem):
        """
        Detect steady state potential
        """
        v_ss = np.max(cell_vmem)
        v_ss_idx = np.argmax(cell_vmem)
        # find_diff = np.diff(cell.vmem)
        return v_ss

    def dV(self, v_ss):
        """ Returns the change in potential dV computed from the steady state
            potential.

        Input: An array containing the steady state potentials from a set of
               simulations

        Returns: 1-D array of change in potential
        """

        # Flattens array of steady state potentials
        v_ss = [i for i in v_ss if i]

        dV = np.zeros(len(v_ss))

        for v in range(len(v_ss)):
            dV[v] = v_ss[v] - self.v_init
        return dV

    def return_dist_to_electrode(self, com_coords):

        self.record_dist = np.zeros(com_coords.shape[0])

        for idx in range(com_coords.shape[0]):
            self.record_dist[idx] = np.sum(np.absolute(
                com_coords[idx, :] - np.array([self.x0, self.y0, self.z0])))

        return self.record_dist

    def run_ext_sim(self, cell, comp_coords, passive=False):
        """ Runs a cell simulation with an added extracellular potential

            Input: LFPy Cell object and coordinates for the comparments to
                   record.
        """

        # print(cell.allsecnames)
        self.create_measure_points(cell, comp_coords)
        # self.print_measure_points(cell)

        self.extracellular_stimuli(cell)
        self.run_cell_simulation(cell)
        # v_max = self.max_mem_pot_dict(cell.vmem)
        self.export_data(cell)

        cell.__del__()
        # return v_max

    def plot_cellsim(self, cell_models_folder, com_coords, morph_ax_params, xlim=[-500, 760], ylim=[-600, 1400], field=False):
        """
        Ploting a combined figure of cell morphology, current amplitude and
        membrane potential
        """

        self.import_data()

        # print(self.cell_tvec[0])
        # Dimensions of the morphology plot
        self.xlim = xlim
        self.ylim = ylim

        # Recreating the cell object used in simulatios without running simul
        cell = self.return_cell(cell_models_folder)
        self.create_measure_points(cell, com_coords)

        # Simulating cell after all parameters and field has been added
        fig = plt.figure(figsize=[10, 8])

        # Defining figure frame and parameters for combined figure
        fig.subplots_adjust(hspace=0.5, left=0.5, wspace=0.5, right=0.96,
                            top=0.9, bottom=0.1)

        self.plot_idxs()

        # Adding morphology to figure
        self.plotSim.plot_morphology(
            cell, fig, xlim, ylim, morph_ax_params)

        n_tsteps = int(self.tstop / self.dt + 1)
        t = np.arange(n_tsteps) * self.dt

        self.insert_pulse(n_tsteps, t)

        if field:
            self.plotSim.plot_external_field(cell, fig)

        # # Setting size and location of plotted potentials and current
        ax_top = 0.90
        ax_h = 0.30
        ax_w = 0.45
        ax_left = 0.5
        stim_axes_placement = [ax_left, ax_top - ax_h, ax_w, ax_h]
        mem_axes_placement = [ax_left, ax_top - ax_h - 0.47, ax_w, ax_h]

        # Adding mambrane potential, current pulse and electrode to combined figure

        self.plotSim.plot_membrane_potential(
            fig, self.cell_tvec, self.cell_vmem, self.tstop, mem_axes_placement)

        self.plotSim.plot_current_pulse(
            fig, self.cell_tvec, self.pulse, self.tstop, stim_axes_placement)

        self.plotSim.draw_electrode(
            self.x0, self.y0, self.z0, self.electrode_radii)

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        self.return_sim_name()
        fig.savefig(join(
            self.save_folder, f'point_source_' + self._sim_name + '.png'))
        # plt.show()
        plt.close(fig=fig)

        cell.__del__()

    def plot_double_morphology(self, cell_models_folder, z_rot, msre_coords, xlim=[-760, 760], ylim=[-500, 1200]):
        """Plots two mrophologies next to each other with an electrode
        representation placed at the same point between them at both morphology
        plots."""

        self.xlim = xlim
        self.ylim = ylim
        morph_ax_params1 = [0.05, 0.1, 0.5, 0.9]
        morph_ax_params2 = [0.473, 0.1, 0.5, 0.9]

        fig = plt.figure()
        self.z_rot = z_rot[0]
        cell1 = self.return_cell(cell_models_folder)
        self.create_measure_points(cell1, msre_coords)
        self.extracellular_stimuli(cell1)

        self.plotSim.plot_idxs(self.measure_pnts)

        self.plotSim.plot_morphology(
            cell1, fig, xlim, ylim, morph_ax_params1)

        if self.x0 > 0:
            self.plotSim.draw_electrode(
                self.x0, self.y0, self.z0, self.electrode_radii)

        cell1.__del__()

        self.z_rot = z_rot[1]
        cell2 = self.return_cell(cell_models_folder)

        self.create_measure_points(cell2, msre_coords)
        # self.extracellular_stimuli(cell2)
        self.plotSim.plot_morphology(
            cell2, fig, xlim, ylim, morph_ax_params2)

        if self.x0 < 0:
            self.plotSim.draw_electrode(
                self.x0, self.y0, self.z0, self.electrode_radii)

        cell2.__del__()

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        fig.savefig(join(
            self.save_folder, f'double_hallermann_{self.cell_name}_z_rot={self.z_rot:.2f}_point_amp={self.amp}uA_elec_pos={self.x0}_{self.y0}_{self.z0}.png'))
        # plt.show()
        plt.close(fig=fig)

        cell2.__del__()

    def plot_double_mem_pot(self):

        self.import_data()

        mem_axes_placement1 = [0.2, 0.6, 0.6, 0.3]
        mem_axes_placement2 = [0.2, 0.2, 0.6, 0.3]
        fig = plt.figure()

        self.plotSim.plot_membrane_potential(
            fig, self.cell_tvec, self.cell_vmem, self.tstop, mem_axes_placement1)

        self.z_rot = np.pi
        self.update_sim_name()
        self.import_data()

        self.plotSim.plot_membrane_potential(
            fig, self.cell_tvec, self.cell_vmem, self.tstop, mem_axes_placement2)

        fig.savefig(join(
            self.save_folder, f'double_vmem_{self.cell_name}_z_rot={self.z_rot:.2f}_point_amp={self.amp}uA_elec_pos={self.x0}_{self.y0}_{self.z0}.png'))
        # plt.show()
        plt.close(fig=fig)
