from Config.Core import Parameters

PDG = Parameters( # named list of PDG identifiers
    down      = Parameters(pdgid = 1),
    up        = Parameters(pdgid = 2),
    strange   = Parameters(pdgid = 3),
    charm     = Parameters(pdgid = 4),
    bottom    = Parameters(pdgid = 5),
    top       = Parameters(pdgid = 6),
    electron  = Parameters(pdgid = 11),
    positron  = Parameters(pdgid = 11),
    muon      = Parameters(pdgid = 13),
    tau       = Parameters(pdgid = 15),
    gluon     = Parameters(pdgid = 21),
    photon    = Parameters(pdgid = 22),
    Z         = Parameters(pdgid = 23),
    W         = Parameters(pdgid = 24),
    proton    = Parameters(pdgid = 2212),
    neutron   = Parameters(pdgid = 2112),
)

def registerParticle(pdgId, name, mass=0., width=0., charge=0, fermion=False):
    PDG[name] = Parameters(
        pdgid = pdgId,
        name = name,
        mass = mass,
        charge = charge,
        width = width,
        fermion = fermion
    )
