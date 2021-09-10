

class NeuronSimulation:
    def __init__():

        self.test = None

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

    def run_cell_simulation(self, cell, vmem=True, imem=False):
        """ Runs the LFPy cell simulation with recordings """
        cell.simulate(rec_vmem=vmem, rec_imem=imem)
