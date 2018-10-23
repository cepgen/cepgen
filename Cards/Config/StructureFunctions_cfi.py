from containers_cfi import Parameters
from SigmaRatio_cfi import SigmaRatio

## Quarks flavours contributions to consider in perturbative SFs definition
class PDFMode:
    ## Consider both valence and sea quarks
    AllQuarks     = 0
    ## Consider only valence quarks
    ValenceQuarks = 1
    ## Consider only sea quarks
    SeaQuarks     = 2

class StructureFunctions:
    '''Types of structure functions supported'''
    generic = Parameters(
        id = 0,
        sigmaRatio = SigmaRatio.E143
    )
    Electron            = generic.clone(id=1)
    ElasticProton       = generic.clone(id=2)
    SuriYennie          = generic.clone(id=11)
    SzczurekUleshchenko = generic.clone(id=12)
    BlockDurandHa       = generic.clone(id=13)
    FioreBrasse         = generic.clone(id=101)
    ChristyBosted       = generic.clone(id=102)
    CLAS                = generic.clone(id=103)
    ALLM91              = generic.clone(id=201)
    ALLM97              = generic.clone(id=202)
    GD07p               = generic.clone(id=203)
    GD11p               = generic.clone(id=204)
    MSTWgrid = generic.clone(
        id = 205,
        gridPath = 'External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat',
    )
    LUXlike = generic.clone(
        id = 301,
        #Q2cut = 10.,
        #W2limits = (4.,1.),
        #continuumSF = GD11p,
        #resonancesSF = ChristyBosted,
    )
    Partonic = generic.clone(
        id = 401,
        pdfSet = 'LUXqed17_plus_PDF4LHC15_nnlo_100',
        numFlavours = 4,
        mode = PDFMode.AllQuarks,
    )
