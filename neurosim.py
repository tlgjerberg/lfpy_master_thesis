from main import ExternalPotentialSimulation
from plotting import PlotSimulation

import numpy as np
import neuron
import LFPy
import os
from os.path import (join, isfile)
from pathlib import Path
import matplotlib.pyplot as plt

font_params = {
    'font.size': 10,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
}

plt.rcParams.update(**font_params)


class NeuronSimulation:
    def __init__(self, cell_params, elec_params):

        self.extPotSim = ExternalPotentialSimulation(elec_params)
        self.plotSim = PlotSimulation()
        self.root_folder = os.path.dirname(__file__)
        self.save_folder = join(
            self.root_folder, cell_params['save_folder_name'])
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
        self.cell_name = cell_params['cell_name']
        self.cell_dist_to_top = cell_params['cell_dist_to_top']
        self.x_shift = cell_params['x_shift']
        self.y_shift = cell_params['y_shift']
        self.y_rot = cell_params['y_rot']
        self.z_rot = cell_params['z_rot']  # Rotation around z-axis

    def return_sim_name(self):
        """Returns a string containing simulation paramters"""

        x = self.extPotSim.x0
        y = self.extPotSim.y0
        z = self.extPotSim.z0
        amp = self.extPotSim.amp
        stop_time = self.extPotSim.stop_time

        self.sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_shift={self.cell_dist_to_top}_z_rot={self.z_rot:.2f}_y_rot={self.y_rot:.2f}_elec_pos={x}_{y}_{z}_t={stop_time}_Amp={amp}mA'

        return self.sim_name

    def return_cell(self, cell_models_folder):
        """Creates a LFPy cell object from a cell model"""

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
                # [ms] dt's should be in powers of 2 for both,
                'dt': self.dt,
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
        # print(name)
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
        """Outputs the compartment name and coordinates of the set of
        record compartments"""

        if hasattr(self, "measure_pnts"):

            print(self.measure_pnts)

            for mc in self.measure_pnts:
                print('name', cell.get_idx_name(mc))
                print('measure_coord', cell.x[mc].mean(),
                      cell.y[mc].mean(), cell.z[mc].mean())

        else:
            raise NameError(
                "create_measure_points() method must be called")

    def run_cell_simulation(self, cell, vmem=True, imem=False):
        """ Runs the LFPy cell simulation with recordings """
        cell.simulate(rec_vmem=vmem, rec_imem=imem)

    def return_tvec(self, cell):

        return cell.tvec

    def return_vmem(self, cell):

        return cell.vmem

    def return_imem(self, cell):

        return cell.imem

    def export_data(self, cell, vmem=True, imem=False):
        """
        Exports arrays of time array, membrane potential and currents of a cells
        simulation to .npy files.

        """

        # Create save folder if it does not exist
        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)
        file_name = self.return_sim_name()

        # Save time array
        np.save(join(self.save_folder, f'{file_name}_tvec'), cell.tvec)

        # Save membrane potential array
        if vmem:
            np.save(join(self.save_folder, f'{file_name}_vmem'), cell.vmem)

        # Save current array
        if imem:
            np.save(join(self.save_folder, f'{file_name}_imem'), cell.imem)

    def import_data(self, vmem=True, imem=False):

        file_name = self.return_sim_name()

        tfile = join(self.save_folder, f'{file_name}_tvec.npy')
        vfile = join(self.save_folder, f'{file_name}_vmem.npy')
        ifile = join(self.save_folder, f'{file_name}_imem.npy')

        if isfile(tfile):
            self.cell_tvec = np.load(tfile)

        else:
            print('File not found!')

        if vmem and isfile(vfile):
            self.cell_vmem = np.load(vfile)

        if imem and isfile(ifile):
            self.cell_imem = np.load(ifile)

    def max_mem_pot_dict(self, cell_vmem):
        """ Returns a dictionary with the maximum membrane potential at each
        recorded compartment chosen as a measurement point. """

        v_max = {}

        for mp in self.measure_pnts:

            v_max[f'{mp}'] = np.max(cell_vmem[mp])

        return v_max

    def run_ext_sim(self, cell, cell_models_folder, comp_coords, passive=False):
        """ Runs a cell simulation with an added extracellular potential """

        # print(cell.allsecnames)
        self.create_measure_points(cell, comp_coords)
        # self.print_measure_points(cell)

        self.extPotSim.extracellular_stimuli(cell, self.tstop, self.dt)
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

    def plot_cellsim(self, cell_models_folder, com_coords, morph_ax_params, xlim=[-500, 760], ylim=[-600, 1400], field=False):
        """
        Ploting a combined figure of cell morphology, current amplitude and
        membrane potential
        """

        self.import_data()

        # Dimensions of the morphology plot
        self.xlim = xlim
        self.ylim = ylim

        # Recreating the cell object used in simulatios without running simul
        cell = self.return_cell(cell_models_folder)
        self.create_measure_points(cell, com_coords)

        self.extPotSim.extracellular_stimuli(cell, self.tstop, self.dt)

        # Simulating cell after all parameters and field has been added
        fig = plt.figure(figsize=[10, 8])

        # Defining figure frame and parameters for combined figure
        fig.subplots_adjust(hspace=0.5, left=0.5, wspace=0.5, right=0.96,
                            top=0.9, bottom=0.1)

        self.plotSim.plot_idxs(self.measure_pnts)

        # Adding morphology to figure
        self.plotSim.plot_morphology(cell, fig, xlim, ylim, morph_ax_params)

        if field:
            self.plot_external_field(cell, fig)

        # # Setting size and location of plotted potentials and current
        ax_top = 0.90
        ax_h = 0.30
        ax_w = 0.45
        ax_left = 0.5
        stim_axes_placement = [ax_left, ax_top - ax_h, ax_w, ax_h]
        mem_axes_placement = [ax_left, ax_top - ax_h - 0.47, ax_w, ax_h]

        # Adding mambrane potential, current pulse and electrode to combined figure
        # self.return_vmem(cell)
        self.plotSim.plot_membrane_potential(
            fig, self.cell_tvec, self.cell_vmem, self.tstop, mem_axes_placement)

        self.plotSim.plot_current_pulse(
            fig, self.cell_tvec, self.extPotSim.pulse, self.tstop, stim_axes_placement)

        x = self.extPotSim.x0
        y = self.extPotSim.y0
        z = self.extPotSim.z0
        electrode_radii = self.extPotSim.electrode_radii
        self.plotSim.draw_electrode(x, y, z, electrode_radii)

        if not os.path.isdir(self.save_folder):
            os.makedirs(self.save_folder)

        self.return_sim_name()
        fig.savefig(join(
            self.save_folder, f'point_source_' + self.sim_name + '.png'))
        # plt.show()
        plt.close(fig=fig)
