from main import ExternalPotentialSim
from parameters import (monophasic_pulse_params, cellsim_bisc_stick_params)
import numpy as np
import neuron
import matplotlib.pyplot as plt
import matplotlib
import LFPy
import os
from os.path import join
import sys
from mpi4py import MPI

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


"""
Improve plotting

Add method for running over a set of electrode positions and current amplitudes

Create plot of Hallermann pyrmidal cell (4 figs?) showing electrode near soma,
dendrite and axon to demonstrate if activation is possible and at what current
amplitude.

Stimulate using current of 10 muA and go lower after (Histed et Al)

Create figure of axial current over distance along axon. Use simple stick model.

Improve plot_cellsim_alt for easy reading and page formatting.

"""
cell_models_folder = join(os.path.dirname(__file__), "cell_models")
current_amps = [-1e4]  # uA
positions = np.array([[0, 0, -50],
                      [0, 0, -100],
                      [0, 0, -200],
                      [0, 0, -400]], dtype=float)
# ,
#[0, 0, -500]
print(positions.shape)
cellsim_bisc_stick_params['save_folder_name'] = 'mpi_axon_test'
axon_measure_idxs = np.array(
    [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])
monophasic_pulse_params['stop_time'] = 200

if RANK == 0:
    sendbuf = positions


else:
    cell_models_folder = None
    current_amps = None  # uA
    positions = None
    axon_measure_idxs = None


COMM.bcast(cell_models_folder, root=0)
COMM.bcast(current_amps, root=0)
COMM.bcast(monophasic_pulse_params, root=0)
COMM.bcast(cellsim_bisc_stick_params, root=0)

recvbuf = np.empty(3, dtype='float')
COMM.Scatter(positions, recvbuf, root=0)

positions = recvbuf
# # Test parameters

# cell_models_folder = join(os.path.dirname(__file__), "cell_models")
# current_amps = [-1e4]  # uA
# positions = np.array([np.array([0, 0, -50], dtype=float),
#                       np.array([0, 0, -100], dtype=float),
#                       np.array([0, 0, -200], dtype=float),
#                       np.array([0, 0, -400], dtype=float),
#                       np.array([0, 0, -500], dtype=float)])

# cellsim_bisc_stick_params['save_folder_name'] = 'mpi_axon_test'
# axon_measure_idxs = np.array(
#     [[0, 0, 0], [0, 0, 300], [0, 0, 600], [0, 0, 1000]])
# monophasic_pulse_params['stop_time'] = 200
print('Rank', RANK)
extPotSim = ExternalPotentialSim(cellsim_bisc_stick_params)

extPotSim.run_ext_sim(cell_models_folder, monophasic_pulse_params, current_amps,
                      positions, axon_measure_idxs, 200, passive=True)
# assert np.allclose(recvbuf, RANK)
