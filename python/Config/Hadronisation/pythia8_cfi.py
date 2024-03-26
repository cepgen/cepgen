##
# \file
# \defgroup pythia8 Pythia 8 hadronisation module
# \ingroup python pythia8
#
# Base Pythia 8 configuration parameters to define the hadronisation module


from EventModifiers.pythia8_cfi import pythia8 as _pythia8
from Config.Hadronisation.pythia8Defaults_cff import pythia8Defaults
from Config.Hadronisation.pythia8Tunes_cff import pythia8CUEP8M1Settings

pythia8 = _pythia8.clone(
    seed = 1000,
    preConfiguration = (
        # printout properties
        # start by disabling some unnecessary output
        'Next:numberCount = 0',
        # disable the decay of resonances by default
        'ProcessLevel:resonanceDecays = off',
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
