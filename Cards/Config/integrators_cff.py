import Config.Core as cepgen
from Config.gsl_cff import GslRngEngine

plain = cepgen.Module('Plain',
    numFunctionCalls = 1000000,
    numPoints = 100,
    rngEngine = GslRngEngine.MT19937,
)

class VegasIntegrationMode:
    Stratified = -1
    ImportanceOnly = 0
    Importance = 1

vegas = plain.clone('Vegas',
    numFunctionCalls = 50000,
    numPoints = 100,
    chiSqCut = 1.5,
    # VEGAS-specific parameters
    iterations = 5,
    alpha = 1.5,
    mode = VegasIntegrationMode.Importance,
    verbosity = -1,
)

miser = plain.clone('MISER',
    # MISER-specific parameters
    estimateFraction = 0.1,
    alpha = 2.,
    dither = 0.,
)
