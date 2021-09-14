import numpy as np
import neuron
import matplotlib
import LFPy
import os
import sys
import time
from os.path import join
from mpi4py import MPI


class ExternalPotentialSimulation:

    def __init__(self, elec_params):

        self.amp = elec_params['pulse_amp']
        self.elec_pos = elec_params['positions']
        self.start_time = elec_params['start_time']
        self.stop_time = elec_params['stop_time']
        self.sigma = elec_params['sigma']
        self.x0, self.y0, self.z0 = elec_params['positions']
        self.pulse_type = elec_params['pulse_type']
        self.electrode_radii = elec_params['electrode_radii']

    def _monophaic_pulse(self, n_tsteps, t):
        # Setting monophasic pulse change at a given time interval
        self.pulse = np.zeros(n_tsteps)
        self.start_idx = np.argmin(np.abs(t - self.start_time))
        self.stop_idx = np.argmin(np.abs(t - self.stop_time))
        self.pulse[self.start_idx:self.stop_idx] = self.amp

    def insert_pulse(self, n_tsteps, t):
        if self.pulse_type == 'monophasic':
            self._monophaic_pulse(n_tsteps, t)

    def extracellular_stimuli(self, cell, tstop, dt):
        """
        Computes the field produced by an extracellular electrode pulse and
        applies it to the cell simulation.

        Parameters: LFPy Cell object, dict of electrode parameters

        Returns:

        Seperate parameter definitions from function?
        """
        # Electrode position

        # Extracellular conductivity

        # Creating external field function
        self.ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi * self.sigma *
                                                           np.sqrt((self.x0 - x)**2 +
                                                                   (self.y0 - y)**2 +
                                                                   (self.z0 - z)**2)))

        # Generating time steps for pulse in simulation
        print(tstop, dt)
        n_tsteps = int(tstop / dt + 1)
        t = np.arange(n_tsteps) * dt

        self.insert_pulse(n_tsteps, t)

        # Applying the external field function to the cell simulation
        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.x.mean(axis=1), cell.y.mean(axis=1),
                                          cell.z.mean(axis=1)).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        # Add the external potential to cell simulation
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
        potential. """

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
