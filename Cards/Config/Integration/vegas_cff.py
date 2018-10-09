from Config.Integration.plain_cff import integrator as plain

class VegasIntegrationMode:
    Stratified = -1
    ImportanceOnly = 0
    Importance = 1

integrator = plain.clone('Vegas',
    numFunctionCalls = 50000,
    chiSqCut = 1.5,
    # VEGAS-specific parameters
    iterations = 10,
    alpha = 1.5,
    mode = VegasIntegrationMode.Importance,
    verbosity = -1,
    loggingOutput = 'cerr',
)

