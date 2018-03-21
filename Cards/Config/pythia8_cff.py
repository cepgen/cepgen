import Config.Core as cepgen

pythia8 = cepgen.Module('pythia8',
    seed = 1000,
    maxTrials = 1,
    pythiaPreConfiguration = (
        # specify we will be using a LHA input
        'Beams:frameType = 5',
        # printout properties
        # start by disabling some unnecessary output
        'Next:numberCount = 0',
        # parameterise the fragmentation part
        #'PartonLevel:Remnants = off',
        # disable all Bremsstrahlung/FSR photon production
        #'PartonLevel:ISR = off',
        #'PartonLevel:FSR = off',
        #'PartonLevel:MPI = off',
        'ParticleDecays:allowPhotonRadiation = off',
        'BeamRemnants:primordialKT = off',
    ),
    pythiaConfiguration = (
        'Tune:preferLHAPDF = 2',
        #'Beams:setProductionScalesFromLHEF = off',
        'SLHA:keepSM = on',
        'SLHA:minMassSM = 1000.',
        'ParticleDecays:limitTau0 = on',
        'ParticleDecays:tau0Max = 10',
        # CUEP8M1 tuning
        'Tune:pp 14', # Monash 2013 tune by Peter Skands (January 2014)
        'MultipartonInteractions:pT0Ref = 2.4024',
        'MultipartonInteractions:ecmPow = 0.25208',
        'MultipartonInteractions:expPow = 1.6',
    ),
    pythiaProcessConfiguration = (),
)
