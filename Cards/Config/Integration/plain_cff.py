from Config.Gsl_cfi import GslRngEngine
from Config.containers_cfi import Module

integrator = Module('plain',
    numFunctionCalls = 1000000,
    rngEngine = GslRngEngine.MT19937,
)
