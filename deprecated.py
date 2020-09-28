def plot_cellsim(self):
    # Simulating cell after all parameters and field has been added
    self.cell.simulate(rec_vmem=True)

    # Setting cell compartments to measure AP
    cell_plot_idxs = [0, int(self.cell.totnsegs / 2),
                      self.cell.totnsegs - 1]
    # cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
    #     1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}
    cell_plot_colors = {idx: ['b', 'cyan', 'orange', 'green',
                              'purple'][num] for num, idx in enumerate(plot_idxs)}

    plt.figure(figsize=(16, 9))

    v_field_ext = np.zeros((50, 200))
    xf = np.linspace(-500, 500, 50)
    zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

    for xidx, x in enumerate(xf):

        for zidx, z in enumerate(zf):
            v_field_ext[xidx, zidx] = self.ext_field(x, 0, z) * self.amp

    vmax = np.max(np.abs(v_field_ext)) / 5
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(121, aspect='equal', xlabel='x [$\mu m$]', ylabel='y [$\mu m$]',
                xlim=[-500, 500], xticks=[-500, 0, 500], title='Green dots: Measurement points')
    plt.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
               origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)

    plt.colorbar(label='mV')
    [plt.plot([self.cell.xstart[idx], self.cell.xend[idx]], [self.cell.zstart[idx], self.cell.zend[idx]], c='gray', zorder=1)
     for idx in range(self.cell.totnsegs)]
    [plt.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o', c=cell_plot_colors[idx], ms=12)
     for idx in cell_plot_idxs]

    l, = plt.plot(self.x0, self.z0, 'y*', ms=2)
    plt.legend([l], ["point current source"], frameon=False)

    # Plotting the membrane potentials
    plt.subplot(222, title='Membrane potential',
                xlabel='Time [ms]', ylabel='mV', ylim=[-80, 20])
    [plt.plot(self.cell.tvec, self.cell.vmem[idx, :], c=cell_plot_colors[idx], lw=2)
     for idx in cell_plot_idxs]

    ax1 = plt.subplot(224, ylim=[-2 * np.max(np.abs(self.pulse / 1000)), 2 * np.max(np.abs(self.pulse / 1000))],
                      ylabel='$\mu$A', title='Injected current')
    ax1.plot(self.cell.tvec, self.pulse / 1000)
    plt.show()


def plot_cellsim(self, measure_idxs):
    cell_plot_idxs = measure_idxs.astype(
        dtype='int')  # List of measurement points
    # cell_plot_colors = {cell_plot_idxs[idx]: plt.cm.Greens_r(
    #     1. / (len(cell_plot_idxs) + 1) * idx + 0.1) for idx in range(len(cell_plot_idxs))}
    cell_plot_colors = idx_clr = {idx: [
        'b', 'cyan', 'orange', 'green', 'purple'][num] for num, idx in enumerate(cell_plot_idxs)}

    # Defining figure frame and parameters
    fig = plt.figure(figsize=[18, 8])
    fig.subplots_adjust(hspace=0.5, left=0.0, wspace=0.5, right=0.96,
                        top=0.9, bottom=0.1)

    # Adding axes with appropriate parameters
    ax_m = fig.add_axes([-0.01, 0.05, 0.2, 0.90], aspect=1, frameon=False,
                        xticks=[], yticks=[], ylim=[-700, 1100], xlim=[-300, 300])

    # Names of different neuron parts and color codings for each
    possible_names = ["Myelin", "axon", "Unmyelin", "Node", "hilloc",
                      "hill", "apic", "dend", "soma"]
    sec_clrs = {"Myelin": 'olive',
                "dend": '0.3',
                "soma": 'k',
                'apic': '0.6',
                "axon": 'lightgreen',
                "Unmyelin": 'salmon',
                "Node": 'r',
                "hilloc": 'lightblue',
                "hill": 'pink', }
    used_clrs = []

    # PLOTTING CELL MORPHOLOGY

    # Sets each segment to the color matching the name set by sec_clrs
    for idx in range(self.cell.totnsegs):
        sec_name = self.cell.get_idx_name(idx)[1]
        # print(sec_name)
        # c = 'k'
        for ax_name in possible_names:
            if ax_name in sec_name:
                # print(ax_name, sec_name)
                c = sec_clrs[ax_name]
                if not ax_name in used_clrs:
                    used_clrs.append(ax_name)

        ax_m.plot([self.cell.xstart[idx], self.cell.xend[idx]],
                  [self.cell.zstart[idx], self.cell.zend[idx]], '-',
                  c=c, clip_on=True, lw=np.sqrt(self.cell.diam[idx]) * 1)

    lines = []
    for name in used_clrs:
        l, = ax_m.plot([0], [0], lw=2, c=sec_clrs[name])
        lines.append(l)
    ax_m.legend(lines, used_clrs, frameon=False,
                fontsize=8, loc=(0.05, 0.0), ncol=2)

    # Plotting dots at the middle of a given section in its given color
    [ax_m.plot(self.cell.xmid[idx], self.cell.zmid[idx], 'o',
               c=cell_plot_colors[idx], ms=13) for idx in cell_plot_idxs]

    ax_m.text(20, 40, "Cortical electrode\n(R={} $\mu$m)".format(self.elec_params["electrode_radii"]),
              fontsize=9, ha='center')

    for e_idx in range(len(self.elec_params["positions"])):
        ellipse_pos = [self.elec_params["positions"][e_idx]
                       [0], self.elec_params["positions"][e_idx][2]]

        ax_m.add_artist(Ellipse(ellipse_pos, width=2 * self.elec_params["electrode_radii"],
                                height=self.elec_params["electrode_radii"] / 5, fc='gray', ec='black'))

    v_field_ext = np.zeros((100, 200))
    xf = np.linspace(-500, 500, 100)
    zf = np.linspace(-500, 1000, 200)
   # print(self.cell.xend)
   # xf = np.linspace(np.min(self.cell.xend), np.max(self.cell.xend), 50)
   # zf = np.linspace(np.min(self.cell.zend), np.max(self.cell.zend), 200)

   for xidx, x in enumerate(xf):

        for zidx, z in enumerate(zf):
            v_field_ext[xidx, zidx] = self.ext_field(
                x, 0, z) * self.amp

    vmax = np.max(np.abs(v_field_ext)) / 5
    ax_cb = plt.gca()
    im_p = ax_cb.imshow(v_field_ext.T, extent=[np.min(xf), np.max(xf), np.min(zf), np.max(zf)],
                        origin='lower', interpolation='nearest', cmap='bwr', vmin=-vmax, vmax=vmax)

    divider = make_axes_locatable(ax_cb)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im_p, cax=cax, label='mV')

    ax_top = 0.90
    ax_h = 0.30
    ax_w = 0.6
    ax_left = 0.3

    ax_vm = fig.add_axes([ax_left, ax_top - ax_h - 0.47, ax_w, ax_h],  # ylim=[-120, 50],
                         xlim=[0, self.tstop], xlabel="Time (ms)")

    ax_vm.set_ylabel("Membrane\npotential (mV)", labelpad=-3)

    # if type(self.spike_time_idxs) == int:
    #     ax_vm.axvline(self.cell.tvec[self.spike_time_idxs], c='r', ls='--')
    ax_stim = fig.add_axes([ax_left, ax_top - ax_h, ax_w, ax_h], xlim=[0, self.tstop],
                           ylabel="Stimuli\ncurrent ($\mu$A)", xlabel="Time (ms)")
    # ax_stim.set_ylabel("$\mu$A", labelpad=-2)
    ax_stim.plot(self.cell.tvec, self.pulse / 1000, lw=0.5)

    # mark_subplots([ax_stim, ax_vm], "BC", xpos=-0.02, ypos=0.98)
    [ax_vm.plot(self.cell.tvec, self.cell.vmem[idx],
                c=cell_plot_colors[idx], lw=0.5) for idx in cell_plot_idxs]

    # plt.show()
    if not os.path.isdir(self.save_folder):
        os.makedirs(self.save_folder)

    fig.savefig(join(
        self.save_folder, f'ext_field_point_amp={self.amp}uA_x={self.x0}_z={self.z0}.png'))
