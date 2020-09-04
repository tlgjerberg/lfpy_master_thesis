

monophasic_pulse_params = dict(pulse_onset=5,
                               pulse_dur=2,
                               # frequency=5000,
                               # asymmetry=-.975,
                               electrode_radii=1,
                               # reverse_pulse_order=True,
                               # pulse_type="biphasic_square",
                               pulse_type="monophasic",
                               # magnitudes=[-0.1e5],
                               pulse_amp=-2.e5,  # Base pulse current amplitude? unit?
                               positions=positions,
                               sigma=0.3,
                               field_type="point_sources")
