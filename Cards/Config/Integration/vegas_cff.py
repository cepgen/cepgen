from Config.Integration.plain_cff import integrator as plain

class VegasIntegrationMode:
    """Type of integration to be performed"""
    Stratified = -1
    ImportanceOnly = 0
    Importance = 1

integrator = plain.clone('Vegas',
    numFunctionCalls = 50000,
    treat = True, # smoothing of the integrand
    chiSqCut = 1.5,
    # VEGAS-specific parameters
    iterations = 10,
    alpha = 1.5,
    mode = VegasIntegrationMode.Importance,
    verbose = -1,
    loggingOutput = 'cerr',
)

