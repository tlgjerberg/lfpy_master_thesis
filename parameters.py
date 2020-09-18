import numpy as np
cortex_top = 0


positions = np.array([[200, 0, -40], ], dtype=float)
cell_dist_to_top = 200
bisc_radius = 8
bisc_elec_pitch = 26  # um
np.random.seed(1234)

monophasic_pulse_params = dict(
    pulse_onset=5,
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
    pulse_amp=1e4,  # nA
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

cellsim_Hallermann_params = dict(
    dt=2**-6,
    tstop=25,
    cut_off=100,
    spike_check_idx_args={'z': -1400, 'section': "my"},
    cell_dist_to_top=cell_dist_to_top,
    z_rot=0,
    x_shift=0,
    cell_name="Hallermann",
    save_folder_name="Hallermann_monophasic",
    save_results=True,
    scale_axon_length=True,
    scale_axon_factor=1.00,
)

cellsim_bisc_Hallermann_params = dict(
    dt=2**-11,
    tstop=20,
    cut_off=20,
    cell_dist_to_top=cell_dist_to_top,
    z_rot=0,
    x_shift=-100,
    cell_name="Hallermann",
    save_folder_name="Hallermann_bisc",
    spike_check_idx_args={'z': -1400, 'section': "my"},
    save_results=True,
    scale_axon_length=True,
    scale_axon_factor=1.25,
)
