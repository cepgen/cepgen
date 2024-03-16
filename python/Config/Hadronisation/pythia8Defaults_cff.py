##
# \file
# \ingroup python pythia8
#
# Collection of Pythia 8 runtime parameters to steer a few decay modes

pythia8Defaults = (
    'Tune:preferLHAPDF = 2',
    #'Beams:setProductionScalesFromLHEF = off',
    #'SLHA:keepSM = on',
    'SLHA:minMassSM = 1000.',
    'ParticleDecays:limitTau0 = on',
    'ParticleDecays:tau0Max = 10',
    'ParticleDecays:allowPhotonRadiation = off',
)
