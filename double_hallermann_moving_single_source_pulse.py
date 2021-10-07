from main import ExternalPotentialSimulation
from plotting import PlotSimulation
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_Hallermann_params)
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
import time
from mpi4py import MPI

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


def run_double_hallermann(cell_models_folder, measure_coords, I, pos, z, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_Hallermann_params['z_rot'] = z

    extPotSim = ExternalPotentialSimulation(
        cellsim_Hallermann_params, monophasic_pulse_params)
    cell = extPotSim.return_cell(cell_models_folder)

    if run_sim:

        extPotSim.run_ext_sim(cell, measure_coords, 20, verbose=True)
        # extPotSim.print_measure_points(cell)

    else:
        print("No simulation run")

    if plot_sim:
        z_rot = [np.pi, 0]
        print('test -1')
        # extPotSim.plot_double_morphology(
        #     cell_models_folder, z_rot, measure_coords)
        print('test 1')
        extPotSim.plot_double_mem_pot(
            z_rot, cell_models_folder, measure_coords)

    else:
        print("No plot generated")


start = time.time()

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_double_morph'

elec_positions = np.array(
    [[700, 60, -352], [750, 60, -352], [-750, -50, -352], [-700, -50, -352]])

I = -1e4

cell_rot = [np.pi, 0]

task_idx = -1

for z in cell_rot:
    if z == 0:
        print(f'{0}')
        measure_coordinates = np.array([[-680, -50, -352]])
    elif z == np.pi:
        print(f'{np.pi:.2f}')
        measure_coordinates = np.array([[680, 60, -352]])
    for pos in elec_positions:
        task_idx += 1
        if not divmod(task_idx, SIZE)[1] == RANK:
            continue

        print('measure coordinates: ', measure_coordinates)
        run_double_hallermann(
            cell_models_folder, measure_coordinates, I, pos, z, False, True)
        print("RANK %d doing task %d" % (RANK, task_idx))

end = time.time()
print(f'Time to execute {end - start} seconds')
