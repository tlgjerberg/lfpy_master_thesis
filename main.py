import numpy as np
import neuron
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import LFPy
import os
import sys
import time
from os.path import join
from matplotlib.patches import Ellipse
from matplotlib.cm import ScalarMappable
from mpi4py import MPI

"""
 transmembrane current i soma/utvalgte punkter som func av tid x
membran pot som func av tid
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

    def return_sim_name(self):
        sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_rot={self.z_rot}_{self.amp}mA_elec_pos={self.elec_pos}'

        return sim_name

    def return_cell(self, cell_models_folder, z=-3, passive=True):

        self.cell_models_folder = cell_models_folder
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
            # self.cell = cell
            self.v_init = cell_parameters['v_init']
            # print(self.cell.zmid[1] - self.cell.zmid[0])
            if not passive:
                neuron.h('forall insert hh')

            return cell

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
            # self.cell = cell
            self.v_init = cell_parameters['v_init']

            # cell.set_pos(z=-np.max(cell.zend) -
            #              self.cell_dist_to_top, x=self.x_shift)
            # cell.set_rotation(z=self.z_rot)
            # self.cell.set_pos(z=-self.cell_dist_to_top)
            # Default rotation
            # self.cell.set_rotation(x=4.729, y=-3.05, z=-3)
            # Axon y-coordinate close to 0
            # self.cell.set_rotation(x=4.9, y=-3.166, z=-3)
            # Apical dendtrite measurement point at aprrox y=0
            # self.cell.set_rotation(x=4.788, y=-3.166, z=-3)
            """
            Adjust positions as needed
            """
            cell.set_pos(z=-self.cell_dist_to_top)
            cell.set_rotation(x=4.729, y=-3.05, z=-3)
            return cell

    def create_measure_points(self, cell, com_coords):
        """
        Setting measurement of membrane potential in compartments closest to
        coordinates.
        """
        measure_pnts = []

        for i in range(com_coords.shape[0]):
            x, y, z = com_coords[i, :]
            measure_pnts.append(cell.get_closest_idx(x, y, z))

        self.measure_pnts = np.array(measure_pnts)

    def _calc_point_sources_field(self, elec_params):
        pass

    def extracellular_stimuli(self, cell, elec_params):
        """
        Parameters: LFPy Cell object, dict of electrode parameters

        Returns:

        Seperate parameter definitions from function?
        """
        self.elec_params = elec_params
        self.elec_pos = elec_params['positions']
        self.amp = elec_params['pulse_amp']  # External current amplitide
        # Electrode position
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
        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.x.mean(axis=1), cell.y.mean(axis=1), cell.z.mean(axis=1)).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        cell.insert_v_ext(v_cell_ext, t)

    def run_cell_simulation(self, cell):
        cell.simulate(rec_vmem=True, rec_imem=True)

    def _find_steady_state(self, cell):
        """
        Detect steady state potential
        """
        self.v_ss = np.max(cell.vmem)
        self.v_ss_idx = np.argmax(cell.vmem)
        find_diff = np.diff(cell.vmem)

    def _dV(self, cell):
        self._find_steady_state(cell)
        self.dV = self.v_ss - self.v_init

    def _record_dist_to_electrode(self, com_coords):

        self.record_dist = np.zeros(com_coords.shape[0])

        for idx in range(com_coords.shape[0]):
            self.record_dist[idx] = np.sum(np.absolute(
                com_coords[idx, :] - np.array([self.x0, self.y0, self.z0])))

        return self.record_dist

    def find_time_constant(self):
        pass

    def run_ext_sim(self, cell_models_folder, elec_params, I,  com_coords, stop_time, elec_pos, idx, passive=False):
        """
        Move cell inside loops of amplitude and postion
        """

        # cell_rot = [-3, -6]
        # for z in cell_rot:
        elec_params['pulse_amp'] = I
        elec_params['positions'] = elec_pos
        cell = self.return_cell(cell_models_folder, passive)
        # self.create_measure_points(cell, com_coords)
        # elec_positions = self.set_electrode_pos(cell, elec_positions)
        elec_dists = np.zeros((len(elec_pos), com_coords.shape[0]))
        ss_pot = np.zeros(len(elec_pos))
        dV = np.zeros(len(elec_pos))

        self.extracellular_stimuli(cell, elec_params)
        self.run_cell_simulation(cell)
        self.export_data(cell)
        # self.plot_cellsim(cell)
        self._find_steady_state(cell)
        ss_pot[idx] = self.v_ss
        self._dV(cell)
        dV[idx] = self.dV
        elec_dists[idx] = self._record_dist_to_electrode(com_coords)

        # self.plot_steady_state(elec_dists[:, 0], ss_pot)
        # self.plot_dV(elec_dists[:, 0], dV)
        # cell.strip_hoc_objects()
        cell.__del__()

    def run_current_sim(self, cell_models_folder, elec_params, current_amps, positions, stop_time, passive=False):
        cell = self.return_cell(cell_models_folder)

        # Neuron activation after cell object has been created
        # if not passive:
        #     neuron.h('forall insert hh')

        elec_params['positions'] = positions
        elec_params['pulse_amp'] = current_amps

        self.extracellular_stimuli(elec_params)
        self.run_cell_simulation()
        # self.plot_currents()

    def export_data(self, cell):
        """
        Export data to text file:
        - compartments measured and their coordinates
        - array of current and membrane potential
        """
        file_name = self.return_sim_name()
        np.save(join(self.save_folder, f'{file_name}_vmem'), cell.vmem)
        np.save(join(self.save_folder, f'{file_name}_tvec'), cell.tvec)
