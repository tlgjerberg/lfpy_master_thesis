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


def run_double_hallermann(cell_models_folder, measure_coords, I, pos, z_rot, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_Hallermann_params['z_rot'] = z_rot[0]

    extPotSim = ExternalPotentialSimulation(
        cellsim_Hallermann_params, monophasic_pulse_params)

    if run_sim:

        cell1 = extPotSim.return_cell(cell_models_folder)
        extPotSim.run_ext_sim(cell1, measure_coords[0], verbose=True)
        cell1.__del__()

        mirror_x_line = 780  # Midpoint
        x_diff = 780 - pos[0]  # Difference between midpoint and electrode
        mr_pos_x = -780 - x_diff

        extPotSim.x0 = mr_pos_x  # Mirrored x coordinate
        extPotSim.y0 = -60
        extPotSim.z_rot = z_rot[1]
        extPotSim.update_sim_name()
        cell2 = extPotSim.return_cell(cell_models_folder)
        extPotSim.run_ext_sim(cell2, measure_coords[1], verbose=True)
        cell2.__del__()
    else:
        print("No simulation run")

    if plot_sim:

        extPotSim.plot_double_morphology(
            cell_models_folder, z_rot, measure_coords)

        extPotSim.plot_double_mem_pot(
            z_rot, cell_models_folder, measure_coords)

    else:
        print("No plot generated")


start = time.time()

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_Hallermann_params['save_folder_name'] = 'data/Hallermann_double_morph'

elec_positions = np.array(
    [[730, 60, -352], [780, 60, -352], [830, 60, -352]])

I = -1e4

cell_rot = [np.pi, 0]

task_idx = -1
measure_coords = np.array([[680, 60, -352], [-680, 60, -352]])

for pos in elec_positions:

    task_idx += 1
    if not divmod(task_idx, SIZE)[1] == RANK:
        continue

    run_double_hallermann(
        cell_models_folder, measure_coords, I, pos, cell_rot, False, True)
    print("RANK %d doing task %d" % (RANK, task_idx))


end = time.time()
print(f'Time to execute {end - start} seconds')
