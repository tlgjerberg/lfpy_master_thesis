from main import ExternalPotentialSimulation
from plotting import PlotSimulation
from fitting import Fitting
from set_electrode_position import set_electrode_pos
from parameters import (monophasic_pulse_params, cellsim_bisc_stick_params)
# from fitting import monoExp, powerlaw, linlaw, fit_exponential, fit_power, fit_linear
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
from mpi4py import MPI
import time
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

"""
Improve plotting

Stimulate using current of 10 muA and go lower after (Histed et Al)

Improve plot_cellsim_alt for easy reading and page formatting.


"""

cell_models_folder = join(os.path.dirname(__file__), "cell_models")
cellsim_bisc_stick_params['save_folder_name'] = 'data/axon_bisc_dist_stim'

# elec_positions = np.array([[0, 0, -50],
#                            [0, 0, -100],
#                            [0, 0, -200],
#                            [0, 0, -400],
#                            [0, 0, -800],
#                            [0, 0, -1600],
#                            [0, 0, -3200]], dtype=float)

# elec_positions = np.array([[0, 0, -200],
#                            [0, 0, -400],
#                            [0, 0, -800],
#                            [0, 0, -1600],
#                            [0, 0, -3200],
#                            [0, 0, -6400],
#                            [0, 0, -12800]], dtype=float)

# elec_positions = np.array([[0, 0, -50],
#                            [0, 0, -100],
#                            [0, 0, -200],
#                            [0, 0, -400],
#                            [0, 0, -800],
#                            [0, 0, -1600],
#                            [0, 0, -3200],
#                            [0, 0, -6400],
#                            [0, 0, -12800]], dtype=float)

elec_positions = np.array([[0, 0, -800],
                           [0, 0, -1600],
                           [0, 0, -3200],
                           [0, 0, -6400],
                           [0, 0, -12800]], dtype=float)


measure_coordinates = np.array(
    [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])

# elec_positions = set_electrode_pos(measure_coordinates)


def run_axon(cell_models_folder, measure_coords, I, pos, z, run_sim=False, plot_sim=False):

    monophasic_pulse_params['pulse_amp'] = I
    monophasic_pulse_params['positions'] = pos
    cellsim_bisc_stick_params['z_rot'] = z

    extPotSim = ExternalPotentialSimulation(
        cellsim_bisc_stick_params, monophasic_pulse_params)
    cell = extPotSim.return_cell(cell_models_folder)

    if run_sim:

        extPotSim.run_ext_sim(cell, measure_coords, z)

        v_ss = extPotSim.find_steady_state_pot(cell.vmem)
    else:
        print('No simulation run!')

    if plot_sim:

        # cell_tvec = np.load(
        #     join(extPotSim.save_folder, f'axon_x_shift=0_z_rot={z:.2f}_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_tvec' + '.npy'))
        # cell_vmem = np.load(
        #     join(extPotSim.save_folder, f'axon_x_shift=0_z_rot={z:.2f}_{I}mA_elec_pos={pos[0]}_{pos[1]}_{pos[2]}_vmem' + '.npy'))
        # plotSim = PlotSimulations(
        #     cellsim_bisc_stick_params, monophasic_pulse_params, cell_vmem, cell_tvec)
        # cell = plotSim.return_cell(cell_models_folder)
        extPotSim.plot_cellsim(
            measure_coords, [0.05, 0.05, 0.3, 0.90], [-300, 300], [-1000, 1100])

    else:
        print('No plots generated!')

    return v_ss


start = time.time()
z = np.pi
I = -1e4  # uA
task_idx = -1
v_ss = None

for pos in elec_positions:
    task_idx += 1
    if not divmod(task_idx, SIZE)[1] == RANK:
        continue

    v_ss = run_axon(cell_models_folder,
                    measure_coordinates, I, pos, z, True, False)

    print("RANK %d doing task %d" % (RANK, task_idx))

# Gathering steady state membrane potentials at root
v_ss = COMM.gather(v_ss, root=0)

if RANK == 0:
    assert v_ss[0] != None
    # Transforming electrode positions into single array of distances
    elec_positions = elec_positions[elec_positions != 0]
    elec_dists = elec_positions.ravel()
    elec_dists *= -1

    extPotSim = ExternalPotentialSimulation(
        cellsim_bisc_stick_params, monophasic_pulse_params)
    cell = extPotSim.return_cell(cell_models_folder)

    dV = extPotSim.dV(v_ss)
    plotSim = PlotSimulation(extPotSim.save_folder)
    plotSim.plot_dV(elec_dists, dV)

    # Fitting the change in potential regular space as a power law
    fit = Fitting(extPotSim.save_folder)
    fit.curve_fit(elec_dists, dV, "powerlaw")
    fit.plot_curve()

    # Fitting the change in potential loglog space as a linear
    fit.curve_fit(np.log(elec_dists), np.log(dV), "linlaw")
    fit.plot_curve()


end = time.time()
print(f'Time to execute {end - start} seconds')
