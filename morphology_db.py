import matplotlib.pyplot as plt
import matplotlib
import os
import zipfile as zip
from os.path import join
from glob import glob
import numpy as np
from neuron import h
import LFPy
from mpi4py import MPI
from unzipper import unzip_directory
import time
import json


COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

with open('layer_download.json') as layer_download:
    layer_data = json.load(layer_download)


def neuron_sets(zip_dir, target_dir, neuron_names):

    if not os.path.isdir(target_dir):
        unzip_directory(zip_dir, target_dir)

    nneurons = []
    for name in neuron_names:
        n = sorted(
            glob(join('hoc_combos_syn.1_0_10.unzipped', 'L5_' + name + '*')))
        nneurons.append(n)
    # Splitting list of neurons into chunks equal to number of processes

    neurons = [j for sub in nneurons for j in sub]

    size = int(len(neurons) / 5)
    neuron_chunks = list(chunks(neurons, size))

    return neuron_chunks


def chunks(l, n):
    """Yield n number of sequential chunks from l."""
    d, r = divmod(len(l), n)
    for i in range(n):
        si = (d + 1) * (i if i < r else r) + d * (0 if i < r else i - r)
        yield l[si:si + (d + 1 if i < r else d)]


Layer_thickness = dict(
    L1=165,
    L23=502,
    L4=190,
    L5=525,
    L6=700
)

Layer_depths = dict(
    L1=165,
    L23=667,
    L4=857,
    L5=1382,
    L6=2082
)


"""
To make sure cell variants are grouped correctly run:
# of cells % # of processes == 5 to to keep the sets of same type of cell in one
batch.
"""
start_time = time.time()

# Create a save folder for data if one does not exist

save_folder = 'data/morphology_search/'

"""
Extract cell models from Blue Brain directory and sort all cells of the same
type into chunks.
"""
if RANK == 0:

    if not os.path.isdir(save_folder):
        os.makedirs(save_folder)

    zip_dir = 'hoc_combos_syn.1_0_10.allzips'
    target_dir = 'hoc_combos_syn.1_0_10.unzipped'

    neuron_names = ['BP', 'BTC', 'DBC', 'MC']  # 6, 3, 7, 7
    # neuron_names = ['STPC', 'TTPC1', 'TTPC2', 'UTPC']  # 1, 1, 1, 1
    neuron_chunks = neuron_sets(zip_dir, target_dir, neuron_names)

else:
    neuron_chunks = None

# Scatter each chunk of neurons between processes
neuron_chunks = COMM.scatter(neuron_chunks, root=0)

"""

"""

save_morph = True
save_cell_hist = True

cell_depths = np.linspace(857, 1382, 100)  # Equally spaced cell depths
cortex_height = 2082  # The total height of the rat cortex in the Blue Brain database


def count_axon_terminals(neuron_chunks, cell_depths, cortex_height):

    terminal_depths = []  # List of axon terminal depths corresonding to indices

    neuron_variants_counter = 0

    # Looping through each set/chunk of neurons and counting the axon terminals
    for nrn in neuron_chunks:
        morphology_path = join(nrn, 'morphology')

        axon_terminals = []  # List of cell indices corresonding to axon terminals

        for morphologyfile in glob(join(morphology_path, '*')):

            morph_name = nrn[31:]  # Name of the morphology

            # Create cell object a given morphology and setting rotation
            cell = LFPy.Cell(morphology=morphologyfile)
            cell.set_rotation(x=(3.0 / 2) * np.pi, y=-0.1)
            for z in cell_depths:
                neuron_variants_counter += 1

                # Loading cell morphology and setting depth
                cell.set_pos(z=z)

                possible_axon_names = ["my", "axon", "ax", "myelin",
                                       "node", "hilloc", "hill"]
                for sec in h.allsec():
                    secname = sec.name()

                    # print(secname)
                    # if h.SectionRef(sec=section).nchild() == 0:

                    # Find index of terminal compartment on the current section
                    extreme_idx = cell.get_idx(section=secname)[-1]

                    # Skip neurons of depths where a neurite falls outside cortex
                    if cell.z[extreme_idx][-1] < 0 or cell.z[extreme_idx][-1] > cortex_height:
                        continue

                    # Test if section is part of an axon
                    sec_is_axon = False
                    for ax_name in possible_axon_names:
                        if ax_name in secname:
                            sec_is_axon = True

                    # We are only interested if sec is axon
                    if not sec_is_axon:
                        continue

                    # If section has no children, then it has an end point
                    if h.SectionRef(sec=sec).nchild() == 0:
                        terminal_idx = cell.get_idx(section=secname)[-1]
                        axon_terminals.append(terminal_idx)
                        terminal_depths.append(
                            (cell.z[terminal_idx].mean(axis=0)))
                        # print('secname: ', secname, 'terminal_idx: ', terminal_idx,
                        #       'z_pos: ', cell.z[terminal_idx].mean(axis=0))

            # Plotting morphology of each cell variant
            if save_morph:
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.invert_yaxis()
                ax1.plot(cell.x.T, cell.z.T, c='k')
                ax1.plot(cell.x[axon_terminals].mean(axis=1),
                         cell.z[axon_terminals].mean(axis=1), 'y*')
                # ax1.set_ylim(cell.z.T[-1], cell.z.T[0])
                ax1.set_xlabel('x[$\mu m$]')
                ax1.set_ylabel('y[$\mu m$]')

                plt.savefig(
                    join(save_folder, f"{morph_name}_morphology_depth.png"), dpi=300)

            return axon_terminals, terminal_depths


axon_terminals, terminal_depths = count_axon_terminals(neuron_chunks)
hist_name = neuron_chunks[0][31:-2]  # Name of neuron type
num_bins = 1000
# Changing terminal depths to shift z-axis to be positive downwards
terminal_depths = np.array(terminal_depths)

# Weights to scale each bin by number of cell variants
const_weights = (1. / len(neuron_chunks)) * np.ones(len(terminal_depths))

# Creating a histogram of axon terminal depths for each neuron type
if save_cell_hist:
    plt.figure()
    plt.hist(terminal_depths, bins=num_bins,
             density=True, orientation='horizontal')
    ax = plt.gca()
    ax.invert_yaxis()
    plt.xlabel('# of terminals')
    plt.ylabel('Layer depth [$\mu m$]')
    plt.savefig(
        join(save_folder, f"layer_terminal_dist_histogram_{hist_name}_bins_{num_bins}.png"), dpi=300)

# Sending length of all terminal_depths arrays to root
sendcounts = np.array(COMM.gather(len(terminal_depths), root=0))

if RANK == 0:
    terminal_depths_layer = np.empty(sum(sendcounts), dtype='float64')

else:
    terminal_depths_layer = None

COMM.Gatherv(terminal_depths, (terminal_depths_layer, sendcounts), root=0)

if RANK == 0:
    print(f'# of terminals in layer {hist_name[0:2]}', len(
        terminal_depths_layer))
    plt.figure()
    # plt.hist(terminal_depths, bins=100,
    #          weights=const_weights, orientation='horizontal')
    plt.hist(terminal_depths_layer, bins=num_bins,
             orientation='horizontal', density=True)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.xlabel('# of terminals')
    plt.ylabel('Layer depth [$\mu m$]')
    plt.savefig(
        join(save_folder, f"layer_terminal_dist_histogram_{hist_name[0:2]}_bins_{num_bins}.png"), dpi=300)

end_time = time.time()
print(f'RANK: {RANK}. Time to execute: {end_time - start_time} seconds')
