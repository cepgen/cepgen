from Config.Integration.plain_cff import integrator as plain

integrator = plain.clone('MISER',
    # MISER-specific parameters
    estimateFraction = 0.1,
    alpha = 2.,
    dither = 0.,
)
