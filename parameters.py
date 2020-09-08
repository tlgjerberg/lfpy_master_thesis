import numpy as np
cortex_top = 0

positions = np.array([[0, 0, cortex_top], ], dtype=float)
cell_dist_to_top = 10
bisc_radius = 8
bisc_elec_pitch = 26  # um
np.random.seed(1234)

monophasic_pulse_params = dict(pulse_onset=5,
                               pulse_dur=2,
                               start_time=2.0,
                               stop_time=5.0,
                               # frequency=5000,
                               # asymmetry=-.975,
                               electrode_radii=1,
                               # reverse_pulse_order=True,
                               # pulse_type="biphasic_square",
                               pulse_type="monophasic",
                               # magnitudes=[-0.1e5],
                               pulse_amp=-2.e5,  # nA
                               positions=positions,
                               sigma=0.3,
                               field_type="point_sources")

cellsim_bisc_stick_params = dict(
    dt=2**-6,
    tstop=30,
    cut_off=10,
    cell_dist_to_top=cell_dist_to_top,  # distance from Stimuli
    # electrode_radii=200,
    z_rot=0,
    x_shift=0,  # Shifting cell in x direction
    cell_name="axon",
    spike_check_idx_args={'z': -550, 'section': "axon"},
    save_folder_name="axon_bisc",
    save_results=True,
)
