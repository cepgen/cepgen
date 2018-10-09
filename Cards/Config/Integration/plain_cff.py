from Config.containers_cfi import Module
from Config.Gsl_cfi import GslRngEngine

integrator = Module('plain',
    numFunctionCalls = 1000000,
    rngEngine = GslRngEngine.MT19937,
)

