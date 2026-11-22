##
# \file
# \ingroup python
# Utilities for the PDG identifiers manipulation

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
# \param pdgid Integer-type PDG identifier
# \param name Computer-safe (and preferably human-readable) particle name
# \param mass Particle on-shell mass (in GeV/c^2)
# \param width Decay width (in GeV)
# \param charge Particle electric charge (in e, may be fractional, e.g. 2./3.)
# \param colour Colour factor (0 for a colourless particle)
# \param fermion Is the particle following the fermion statistics?
# \note The C++ particle properties store the electric charge as a list of integers, in units of e/3, under the key
#  'charges', and the colour factor under the key 'colours'.
def registerParticle(pdgid: int,
                     name: str,
                     mass: float=0.,
                     width: float=0.,
                     charge: float=0.,
                     colour: int=1,
                     fermion: bool=False):
    globals()['PDG'][name] = Parameters(
        name = name,
        description = name,
        pdgid = pdgid,
        mass = mass,
        charges = (int(round(3. * charge)),),
        width = width,
        colours = colour,
        fermion = fermion
    )
    print('particle with pdg={} defined: {}'.format(pdgid, globals()['PDG'][name]))
