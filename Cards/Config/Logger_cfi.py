from .containers_cfi import Parameters

class Logging:
    '''Logging verbosity'''
    Nothing         = 0
    Error           = 1
    Warning         = 2
    Information     = 3
    Debug           = 4
    DebugInsideLoop = 5

logger = Parameters(
    level = Logging.Information,
    enabledModules = ('',),
)
