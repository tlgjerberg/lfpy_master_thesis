import numpy as np
import neuron
import matplotlib
import LFPy
import os
import sys
import time
from os.path import join
from mpi4py import MPI

"""
 transmembrane current i soma/utvalgte punkter som func av tid x
membran pot som func av tid
"""


class ExternalPotentialSim:

    def __init__(self, cell_params, elec_params):

        self.root_folder = os.path.dirname(__file__)
        self.cell_params = cell_params
        self.elec_params = elec_params
        self.save_folder = join(
            self.root_folder, cell_params['save_folder_name'])
        self.amp = elec_params['pulse_amp']
        self.elec_pos = elec_params['positions']
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
        self.cell_name = cell_params['cell_name']
        self.cell_dist_to_top = cell_params['cell_dist_to_top']
        self.x_shift = cell_params['x_shift']
        self.y_shift = cell_params['y_shift']
        self.y_rot = cell_params['y_rot']
        self.z_rot = cell_params['z_rot']  # Rotation around z-axis
        self.start_time = elec_params['start_time']
        self.stop_time = elec_params['stop_time']

    def return_sim_name(self):
        self.sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_shift={self.cell_dist_to_top}_z_rot={self.z_rot:.2f}_y_rot={self.y_rot:.2f}_elec_pos={self.elec_pos[0]}_{self.elec_pos[1]}_{self.elec_pos[2]}_t={self.stop_time}_Amp={self.amp}mA'

        return self.sim_name

    def return_cell(self, cell_models_folder):

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
            self.v_init = cell_parameters['v_init']
            cell = LFPy.Cell(**cell_parameters)
            cell.set_rotation(x=0, y=self.y_rot, z=self.z_rot)
            cell.set_pos(x=self.x_shift, y=self.y_shift,
                         z=self.cell_dist_to_top)
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

            self.v_init = cell_parameters['v_init']
            cell = LFPy.Cell(**cell_parameters)

            # Adjusting the cell position and rotation
            cell.set_pos(x=self.x_shift, y=self.y_shift,
                         z=-self.cell_dist_to_top)
            cell.set_rotation(x=4.729, y=-3.05, z=self.z_rot)
            return cell

    def find_terminals(self, cell, measure_pnt):

        idx = cell.get_closest_idx(measure_pnt)
        name = cell.get_idx_name(idx=idx)
        print(name)
        # term_idx = cell.get_idx_children(parent=)
        #
        # if h.SectionRef(sec=sec).nchild() == 0:
        #     pass

    def create_measure_points(self, cell, com_coords):
        """
        Setting measurement of membrane potential in compartments closest to
        coordinates.
        """

        measure_pnts = []

        if com_coords.ndim == 1:

            x, y, z = com_coords
            measure_pnts.append(cell.get_closest_idx(x, y, z))

        else:

            for i in range(com_coords.shape[0]):

                x, y, z = com_coords[i, :]
                measure_pnt = cell.get_closest_idx(x, y, z)

                self.find_terminals(cell, measure_pnt)
                measure_pnts.append(measure_pnt)

        self.measure_pnts = np.array(measure_pnts)

    def print_measure_points(self, cell):

        if hasattr(self, "measure_pnts"):

            print(self.measure_pnts)

            for mc in self.measure_pnts:
                print('name', cell.get_idx_name(mc))
                print('measure_coord', cell.x[mc].mean(),
                      cell.y[mc].mean(), cell.z[mc].mean())

        else:
            raise NameError("create_measure_points() method must be called")

    def extracellular_stimuli(self, cell):
        """
        Parameters: LFPy Cell object, dict of electrode parameters

        Returns:

        Seperate parameter definitions from function?
        """
        # Electrode position
        # self.elec_pos = self.elec_params['positions']
        self.x0, self.y0, self.z0 = self.elec_params['positions']

        sigma = self.elec_params['sigma']

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
        self.pulse[self.start_idx:self.stop_idx] = self.amp

        # Applying the external field function to the cell simulation
        v_cell_ext = np.zeros((cell.totnsegs, n_tsteps))
        v_cell_ext[:, :] = self.ext_field(cell.x.mean(axis=1), cell.y.mean(axis=1),
                                          cell.z.mean(axis=1)).reshape(
            cell.totnsegs, 1) * self.pulse.reshape(1, n_tsteps)

        cell.insert_v_ext(v_cell_ext, t)

    def run_cell_simulation(self, cell, vmem=True, imem=False):
        cell.simulate(rec_vmem=vmem, rec_imem=imem)

    def find_steady_state_pot(self, cell_vmem):
        """
        Detect steady state potential
        """
        v_ss = np.max(cell_vmem)
        v_ss_idx = np.argmax(cell_vmem)
        # find_diff = np.diff(cell.vmem)
        return v_ss

    # def find_max_mem_pot(self, cell_vmem, measure_pnts_idx):
    #
    #
    #     for mp in self.measure_pnts:
    #
    #
    #         v_max = np.max(cell_vmem[mp])
    #
    #     return v_max

    def max_mem_pot_dict(self, cell_vmem):

        v_max = {}

        for mp in self.measure_pnts:

            v_max[f'{mp}'] = np.max(cell_vmem[mp])

        return v_max

    def dV(self, v_ss):
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

    def return_segment_coords(self, cell):

        for idx in self.measure_pnts:
            coords = [cell.x[idx].mean(axis=0), cell.y[idx].mean(
                axis=0), cell.z[idx].mean(axis=0)]

    def export_data(self, cell, vmem=True, imem=False):
        """
        Export data to text file:
        - compartments measured and their coordinates
        - array of current and membrane potential
        """
        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)
        file_name = self.return_sim_name()

        np.save(join(self.save_folder,
                     f'{file_name}_tvec'), cell.tvec)

        np.save(join(self.save_folder,
                     f'{file_name}_vmem'), cell.vmem)

        if imem:
            np.save(join(self.save_folder,
                         f'{file_name}_imem'), cell.imem)

    def run_ext_sim(self, cell_models_folder, comp_coords, passive=False):

        cell = self.return_cell(cell_models_folder)
        # print(cell.allsecnames)
        self.create_measure_points(cell, comp_coords)
        self.print_measure_points(cell)
        self.return_segment_coords(cell)

        self.extracellular_stimuli(cell)
        self.run_cell_simulation(cell)
        v_max = self.max_mem_pot_dict(cell.vmem)
        self.export_data(cell)

        cell.__del__()
        return v_max

    def run_current_sim(self, cell_models_folder, comp_coords, passive=False):
        cell = self.return_cell(cell_models_folder)
        self.create_measure_points(cell, comp_coords)
        self.extracellular_stimuli(cell)
        self.run_cell_simulation(cell, False, True)
        self.export_data(cell, False, True)
