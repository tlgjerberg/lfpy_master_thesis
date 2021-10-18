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
from collections import defaultdict
from random import randint

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


def neuron_sets(zip_dir, target_dir, neuron_names):
    """
    Extract cell models from Blue Brain directory and sort all cells of the same
    type into chunks.
    """
    if not os.path.isdir(target_dir):

        unzip_directory(zip_dir, target_dir)

    neuron_chunks = []
    neuron_chunks_dict = {}
    for name in neuron_names:
        neuron_chunks_dict = {}

        n = sorted(
            glob(join('hoc_combos_syn.1_0_10.unzipped', name + '*')))

        neuron_chunks_dict[name] = n
        neuron_chunks.append(neuron_chunks_dict)

    return neuron_chunks


def count_axon_terminals(neuron_type, neuron_sub_types, cell_depths, cortex_height, save_morph=False):
    """
    To make sure cell variants are grouped correctly run:
    # of cells % # of processes == 5 to to keep the sets of same type of cell in one
    batch.
    """

    terminal_depths = []  # List of axon terminal depths corresonding to indices

    neuron_variants_counter = 0

    # Looping through each set/chunk of neurons and counting the axon terminals
    for nrn in neuron_sub_types:
        morphology_path = join(nrn, 'morphology')

        axon_terminals = []  # List of cell indices corresonding to axon terminals

        morph_name = nrn[31:-2]  # Name of the morphology

        for morphologyfile in glob(join(morphology_path, '*')):

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


Layer_thickness = dict(
    L1=165,
    L23=502,
    L4=190,
    L5=525,
    L6=700
)

Layer_depths = dict(
    L1=0,
    L23=165,
    L4=667,
    L5=857,
    L6=1382
)


start_time = time.time()

# Create a save folder for data if one does not exist
save_folder = 'data/morphology_search/'


with open('layer_download.json') as layer_download:
    layer_data = json.load(layer_download)

layer_L5 = layer_data['L5']  # Dictionary of layer 5 data
# Layer 5 morphological type data
L5_morph_types = layer_L5['No. of neurons per morphological types']

if RANK == 0:

    if not os.path.isdir(save_folder):
        os.makedirs(save_folder)

    zip_dir = 'hoc_combos_syn.1_0_10.allzips'
    target_dir = 'hoc_combos_syn.1_0_10.unzipped'

    # neuron_names = ['L5_BP', 'L5_BTC', 'L5_DBC', 'L5_MC']  # 6, 3, 7, 7
    neuron_names = ['L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC']  # 1, 1, 1, 1
    neuron_chunks = neuron_sets(zip_dir, target_dir, neuron_names)

else:
    neuron_chunks = None

# Scatter each chunk of neurons between processes
neuron_chunks = COMM.scatter(neuron_chunks, root=0)


cell_depths = np.linspace(857, 1382, 100)  # Equally spaced cell depths in L5
# The total height of the rat cortex in the Blue Brain database
cortex_height = Layer_depths['L6'] + Layer_thickness['L6']

# Split dictionary of neuron type and morphology files
neuron_type, neuron_sub_types = next(
    iter(neuron_chunks.items()))

axon_terminals, terminal_depths = count_axon_terminals(neuron_type,
                                                       neuron_sub_types,
                                                       cell_depths, cortex_height)

# Creating a histogram of axon terminal depths for each neuron subtype
subtype_name = neuron_sub_types[0][31:-2]  # Subtype name

num_bins = 1000  # Number of bins used in the histograms

plt.figure()
plt.hist(terminal_depths, bins=num_bins,
         density=True, orientation='horizontal')
ax = plt.gca()
ax.invert_yaxis()
plt.xlabel('# of terminals')
plt.ylabel('Layer depth [$\mu m$]')
plt.savefig(
    join(save_folder, f"layer_terminal_dist_histogram_{subtype_name}_bins_{num_bins}.png"), dpi=300)


terminal_depths_dict_list = {neuron_type: terminal_depths}
print(neuron_type, len(terminal_depths))
terminal_depths_dict_list = COMM.gather(terminal_depths_dict_list, root=0)


if RANK == 0:

    terminal_depths_dict = defaultdict(list)
    for d in terminal_depths_dict_list:
        for key, value in d.items():
            terminal_depths_dict[key].extend(value)

    print(type(terminal_depths_dict))

    num_bins = 1000  # Number of bins used in the histograms

    for key, value in terminal_depths_dict.items():
        print(key, len(value))

        # Number of neurons of selected type
        num_neurons = L5_morph_types[neuron_type]
        scale_factor = 1. / num_neurons  # Weights scaling
        type_weights = scale_factor * \
            np.ones(len(terminal_depths_dict[key]))  # Array of weights

        # Creating a histogram of axon terminal depths for each neuron type
        plt.figure()
        plt.hist(value, bins=num_bins,
                 density=True, weights=type_weights, orientation='horizontal')
        ax = plt.gca()
        ax.invert_yaxis()
        plt.xlabel('# of terminals')
        plt.ylabel('Layer depth [$\mu m$]')
        plt.savefig(
            join(save_folder, f"layer_terminal_dist_histogram_{key}_bins_{num_bins}.png"), dpi=300)

    # plt.figure()
    # plt.hist(terminal_depths_layer, bins=num_bins,
    #          orientation='horizontal', density=True)
    # ax = plt.gca()
    # ax.invert_yaxis()
    # plt.xlabel('# of terminals')
    # plt.ylabel('Layer depth [$\mu m$]')
    # plt.savefig(
    #     join(save_folder, f"layer_terminal_dist_histogram_{neuron_type}_bins_{num_bins}.png"), dpi=300)

end_time = time.time()
print(f'RANK: {RANK}. Time to execute: {end_time - start_time} seconds')
