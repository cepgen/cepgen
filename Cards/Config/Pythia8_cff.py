from containers_cfi import Module, Parameters
from pythia8Defaults_cfi import pythia8Defaults
from pythia8Tunes_cfi import pythia8CUEP8M1Settings

pythia8 = Module('pythia8',
    moduleParameters = Parameters(
        seed = 1000,
        maxTrials = 1,
    ),
    preConfiguration = (
        # printout properties
        # start by disabling some unnecessary output
        'Next:numberCount = 0',
        # parameterise the fragmentation part
        #'PartonLevel:Remnants = off',
        # disable all Bremsstrahlung/FSR photon production
        'PartonLevel:ISR = off',
        'PartonLevel:FSR = off',
        'PartonLevel:MPI = off',
        'BeamRemnants:primordialKT = off',
    ),
    commonSettings = pythia8Defaults,
    tuningSettings = pythia8CUEP8M1Settings,
    #tuningSettings = pythiaCUEP8S1CTEQ6L1Settings,
    #tuningSettings = pythiaCUEP8M2T4Settings,
    processConfiguration = (
        'commonSettings', 'tuningSettings',
    ),
)
