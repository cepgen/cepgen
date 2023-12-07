from .SigmaRatio_cfi import SigmaRatio
from .containers_cfi import Module


## Quarks flavours contributions to consider in perturbative SFs definition
class PDFMode:
    ## Consider both valence and sea quarks
    AllQuarks     = 0
    ## Consider only valence quarks
    ValenceQuarks = 1
    ## Consider only sea quarks
    SeaQuarks     = 2

class StructureFunctions:
    """Types of structure functions supported"""
    generic = Module(0,
        sigmaRatio = SigmaRatio.SibirtsevBlunden
    )
    Electron            = generic.clone(1)
    ElasticProton       = generic.clone(2)
    SuriYennie          = generic.clone(11)
    SzczurekUleshchenko = generic.clone(12)
    BlockDurandHa       = generic.clone(13)
    SuriYennieAlt       = generic.clone(14)
    FioreBrasse         = generic.clone(101)
    ChristyBosted       = generic.clone(102)
    CLAS                = generic.clone(103)
    FioreBrasseLowMass  = generic.clone(104)
    DonnachieLandshoff  = generic.clone(105)
    Capella             = generic.clone(106)
    ALLM91              = generic.clone(201)
    ALLM97              = generic.clone(202)
    GD07p               = generic.clone(203)
    GD11p               = generic.clone(204)
    MSTWgrid = generic.clone(205,
        gridPath = 'External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat',
    )
    HHT                 = generic.clone(206)
    HHT_FT              = generic.clone(207)
    LUXlike = generic.clone(301,
        #Q2cut = 10.,
        #W2limits = (4.,1.),
        #continuumSF = GD11p,
        #resonancesSF = ChristyBosted,
    )
    Shamov              = generic.clone(302)
    KulaginBarinov      = generic.clone(303)
    BodekKangXu         = generic.clone(304)
    Partonic = generic.clone(401,
        pdfSet = 'LUXqed17_plus_PDF4LHC15_nnlo_100',
        numFlavours = 4,
        mode = PDFMode.AllQuarks,
    )
    APFELPartonic       = generic.clone(402)
