from Config.containers_cfi import Module, Parameters

class DecayMode:
    All = 0
    ElectronMode = 1
    MuonMode = 2
    PionMode = 3
    RhoMode = 4
    A1Mode = 5
    KMode = 6
    KStarMode = 7
    FourPionMode = 8 # 2 pi^pm pi^mp pi^0 nu
    FourPionNeutralMode = 9 # 3 pi^0 pi^pm nu
    FivePionNeutralMode = 10 # 2 pi^pm pi^mp 2 pi^0 nu
    FivePionMode = 11 # 3 pi^pm 2 pi^mp nu
    SixPionMode = 12 # 3 pi^pm 2 pi^mp pi^0 nu
    SixPionNeutralMode = 13 # 2 pi^pm pi^mp 3 pi^0 nu
    TwoKPionMode = 14 # K^pm K^mp pi^pm nu
    TwoKPionNeutralMode = 15 # K^0 Kbar^0 pi^pm nu
    TwoKNeutralPionMode = 16 # K^pm K^0 pi^0 nu
    TwoPionNeutralKMode = 17 # 2 pi^0 K^pm nu
    TwoPionKMode = 18 # pi^pm pi^mp K^pm nu
    TwoPionKNeutralMode = 19 # pi^pm pi^0 Kbar^0 nu
    TwoPionEtaMode = 20 # eta pi^pm pi^0 nu
    TwoPionGammaMode = 21 # pi^pm pi^0 gamma nu
    TwoKMode = 22 # K^pm K^0 nu

tauola = Module('tauola',
    polarisations = Parameters(
        full = True,
        #GAMMA = False,
    ),
)
