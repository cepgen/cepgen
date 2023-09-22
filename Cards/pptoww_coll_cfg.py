import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cff import generator as _gen

process = coll.process.clone('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        method = 0,  # on-shell (0) or off-shell (1) formula
        polarisationStates = 0,  # full
    ),
    inKinematics = cepgen.Parameters(
        #partonFluxes = (coll.ProtonFlux.LHAPDFLUXlep, coll.ProtonFlux.LHAPDFLUXlep),
        partonFluxes = (coll.ProtonFlux.PhotonElastic, coll.ProtonFlux.PhotonElastic),
        cmEnergy = 13.e3,
    ),
    outKinematics = coll.process.outKinematics.clone(
        mx = (1.07, 1000.),
        q2 = (0., 10.),
        #--- extra cuts on the W+W- system
        invmass = (0.,),
        #--- cuts on single particles' level
        cuts = {
            # cuts on the single W level
            #24: cepgen.Parameters(pt = (0.,)), # no pt cut on Ws
            # cuts on the W decay products
            # (mimicking LHC-like experimental cuts)
            #11: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5)),
            #13: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5))
        },
        #xi = (0.02, 0.15),
    )
)

generator = _gen.clone(
    numEvents = 50000,
    printEvery = 5000,
)
text = cepgen.Module('text',  # histogramming/ASCII output capability
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(50., 1000.), nbins=19),
        'm(ob2)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'acop(7,8)': cepgen.Parameters(nbins=10, log=True),
    }
)
output = cepgen.Sequence(text)
