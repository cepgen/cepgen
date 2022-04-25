"""\file
\brief Utilities for the PDG identifiers manipulation
"""

from Config.Core import Parameters

## Named list of PDG identifiers
# \note Minimal information required is the pdgid (unsigned int) attribute
PDG = Parameters(
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

## Define a new particle type into the PDG library
# \param pdgId Integer-type PDG identifier
# \param name Computer-safe (and preferably human-readable) particle name
# \param mass Particle on-shell mass (in GeV/c^2)
# \param width Decay width (in GeV)
# \param charge Particle electric charge (in e)
# \param colour Colour charge
# \param fermion Is the particle following the fermion statistics?
def registerParticle(pdgId, name, mass=0., width=0., charge=0, colour=1, fermion=False):
    PDG[name] = Parameters(
        name = name,
        description = name,
        pdgid = pdgId,
        mass = mass,
        charge = charge,
        width = width,
        colour = colour,
        fermion = fermion
    )
