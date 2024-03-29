from plotting import PlotSimulation
from glob import glob
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
    """
    A class for simulating a neuron using LFPy and NEURON and extracting features
    and inforamtion.

    Attributes:
    ----------
    root_folder : str
        Path to working directory
    save_folder : str
        Path to save directory
    dt : float
        Simulation time step
    tstop : int
        Simulation end
    cut_off : float
        Simulation start time
    cell_name : str
        Name of cell model
    cell_dist_to_top : int
        z-direction coordinate shift
    x_shift : int
        x-direction coordinate shift
    y_shift : int
        y-direction coordinate shift
    y_rot : float
        Rotation around y-axis
    z_rot : float
        Rotation around y-axis

    Methods:
    -------
    return_sim_name()

    return_cell(cell_models_folder)

    find_terminals(cell, measure_pnt)

    create_measure_points(cell, com_coords)

    find_secnames(self, cell)

    print_measure_points(cell)

    run_cell_simulation(cell, vmem=True, imem=False)

    export_data(cell, vmem=True, imem=False)

    import_data(vmem=True, imem=False)

    max_mem_pot_dict(cell_vmem)

    consolidate_v_max(sim_name, idx)

    run_current_sim(cell_models_folder, comp_coords, passive=False)
    """

    def __init__(self, cell_params):

        self.root_folder = os.path.dirname(__file__)
        self.save_folder = join(
            self.root_folder, cell_params['save_folder_name'])
        self.dt = cell_params['dt']
        self.tstop = cell_params['tstop']
        self.cut_off = cell_params["cut_off"]
        self.cell_name = cell_params['cell_name']  # Name of cell model
        self.cell_dist_to_top = cell_params['cell_dist_to_top']
        self.x_shift = cell_params['x_shift']  # x-direction coordinate shift
        self.y_shift = cell_params['y_shift']  # y-direction coordinate shift
        self.y_rot = cell_params['y_rot']  # Rotation around y-axis
        self.z_rot = cell_params['z_rot']  # Rotation around z-axis

        self.plotSim = PlotSimulation(self.save_folder)
        self._sim_name = f'{self.cell_name}_x_shift={self.x_shift}_z_shift={self.cell_dist_to_top}_z_rot={self.z_rot:.2f}_y_rot={self.y_rot:.2f}'

    def return_sim_name(self):
        """Returns a string containing simulation parameters"""

        return self._sim_name

    def return_cell(self, cell_models_folder):
        """
        Creates a LFPy cell object from a cell model

        Parameters:
        cell_models_folder (str): Name of the directory location of cell models

        Returns:
        cell: LFPy Cell object determined by model and parameters
        """

        if self.cell_name == 'axon':

            model_path = join(cell_models_folder, 'unmyelinated_axon.hoc')

            cell_parameters = {
                'morphology': model_path,
                'nsegs_method': "lambda_f",
                'lambda_f': 1000.,
                'v_init': -65,  # initial crossmembrane potential
                'passive': False,  # Toggle passive mechanisms
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
                         z=self.cell_dist_to_top - cell.z[0][0])

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
                'passive': False,   # Toggle passive mechanisms
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

        Parameters:
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

    def find_secnames(self, cell):
        """
        Finds the name of a section matching each measurement compartment and
        adds both to a dictionary as a value and key pair.

        Input:
        cell: LFPy Cell object
        """

        self.sec_names = {}

        for mp in self.measure_pnts:

            idx, secname, _ = cell.get_idx_name(mp)

            self.sec_names[f'{idx}'] = secname

    def print_measure_points(self, cell):
        """Outputs the compartment name and coordinates of the set of
        record compartments"""

        if hasattr(self, "measure_pnts"):

            print(self.measure_pnts)

            for mc in self.measure_pnts:
                print('name', cell.get_idx_name(mc))
                print('measurement coordinate', cell.x[mc].mean(),
                      cell.y[mc].mean(), cell.z[mc].mean())

        else:
            raise NameError(
                "create_measure_points() method must be called")

    def run_cell_simulation(self, cell, vmem=True, imem=False):
        """ Runs the LFPy cell simulation with recordings """
        cell.simulate(rec_vmem=vmem, rec_imem=imem)

    @property
    def return_tvec(self, cell):

        return cell.tvec

    @property
    def return_vmem(self, cell):

        return cell.vmem

    @property
    def return_imem(self, cell):

        return cell.imem

    def export_data(self, cell, vmem=True, imem=False):
        """
        Exports arrays of time array, membrane potential and currents of a cells
        simulation to .npy files.

        Parameters:
        cell (LFPy.Cell): LFPy Cell object
        vmem (bool): Arguement to export membrane potential file
        imem (bool): Arguement to export transmembrane current file
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
        """
        Imports data from simulations saved to file

        Parameters:
        vmem (bool): Arguement to import membrane potential file
        imem (bool): Arguement to import transmembrane current file
        """

        file_name = self._sim_name

        tfile = join(self.save_folder, f'{file_name}_tvec.npy')
        vfile = join(self.save_folder, f'{file_name}_vmem.npy')
        ifile = join(self.save_folder, f'{file_name}_imem.npy')

        if isfile(tfile):
            self.cell_tvec = np.load(tfile)

        else:
            print('File not found. Run Simulation!')
            raise FileNotFoundError

        if vmem and isfile(vfile):
            self.cell_vmem = np.load(vfile)

        if imem and isfile(ifile):
            self.cell_imem = np.load(ifile)

    def max_mem_pot_dict(self, cell_vmem):
        """
        Returns a dictionary with the maximum membrane potential at each
        recorded compartment chosen as a measurement point.

        Parameters:
        cell_vmem (array): Cell membrane potential recordings

        Returns:
        v_max: Dictionary of maximum potential for each measured compartment
        """
        v_max = {}

        for mp in self.measure_pnts:

            v_max[f'{mp}'] = np.max(cell_vmem[mp])

        return v_max

    def consolidate_v_max(self, sim_name, idx):
        """
        Loads a set of memembrane potential simulations and consolidates a list of
        all maximum potentials.

        Parameters:
        sim_name (str): Path to a set of simulations.
        idx (int): index of the compartment of interest.

        Returns:
        vml: A list of maximum potentials.
        """

        vmem_list = sorted(glob(
            join(self.save_folder, sim_name + '*_vmem.npy')))

        if vmem_list:
            print('Files found')

        else:
            raise FileNotFoundError

        vml = []
        for vmem in vmem_list:
            cell_vmem = np.load(vmem)
            v_max = np.max(cell_vmem[idx])
            vml.append(v_max)

        return vml

    def run_current_sim(self, cell_models_folder, comp_coords, passive=False):
        """
        Runs a set of methods for simulating the cell and and marking a spesific
        set of compartments

        Parameters:
        cell (LFPy.Cell): LFPy Cell object determined by model and parameters

        """
        cell = self.return_cell(cell_models_folder)
        self.create_measure_points(cell, comp_coords)
        self.extracellular_stimuli(cell)
        self.run_cell_simulation(cell, False, True)
        self.export_data(cell, False, True)
