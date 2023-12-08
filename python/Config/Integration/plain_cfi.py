from Config.Gsl_cff import GslRngEngine
from Config.containers_cff import Module

plain = Module('plain',
    numFunctionCalls = 1000000,
    rngEngine = GslRngEngine.MT19937,
)
