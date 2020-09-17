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
