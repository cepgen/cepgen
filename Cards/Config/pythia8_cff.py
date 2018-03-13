import Config.Core as cepgen

pythia8 = cepgen.Module('pythia8',
    seed = 0,
    maxTrials = 1,
    pythiaPreConfiguration = (
        #'Init:showAllSettings = on',
        # disable all generation processes
        #'ProcessLevel:all = off',
        # printout properties
        # start by disabling some unnecessary output
        'Next:numberCount = 0',
        #'Next:numberShowInfo = 0',
        #'Next:numberShowProcess = 0',
        # parameterise the fragmentation part
        #'BeamRemnants:beamJunction = on',
        #'ProcessLevel:resonanceDecays = on',
        #'PartonLevel:Remnants = on',
        #'PartonLevel:all = off',
        #'PartonLevel:MPI = on',
        #'HadronLevel:all = on',
        #'Diffraction:doHard = on',
        #'HardQCD:all = on',
        # disable all Bremsstrahlung/FSR photon production
        #'PartonLevel:all = on',
        'PartonLevel:ISR = off',
        'PartonLevel:FSR = off',
        'PartonLevel:MPI = on',
        'ParticleDecays:allowPhotonRadiation = off',
    ),
    pythiaConfiguration = (
        'Tune:preferLHAPDF = 2',
        #'Main:timesAllowErrors = 10000',
        #'Check:epTolErr = 0.01',
        'Check:abortIfVeto = off',
        'Beams:setProductionScalesFromLHEF = off',
        'SLHA:keepSM = on',
        'SLHA:minMassSM = 1000.',
        'ParticleDecays:limitTau0 = on',
        'ParticleDecays:tau0Max = 10',
        # CUEP8M1 tuning
        'Tune:pp 14', # Monash 2013 tune by Peter Skands (January 2014)
        'MultipartonInteractions:pT0Ref = 2.4024',
        'MultipartonInteractions:ecmPow = 0.25208',
        'MultipartonInteractions:expPow = 1.6',
        #'Check:epTolErr = 10000.', #FIXME FIXME FIXME!!! BAD BAD BAD
        #'Check:epTolWarn = 10.', #FIXME FIXME FIXME!!! BAD BAD BAD
        '24:onMode = off', # disable all W decays...
        '13:onMode = off', # disable all muon decays...
    ),
    pythiaProcessConfiguration = (),
)
