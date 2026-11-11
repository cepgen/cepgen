##
# \file
# \defgroup sherpa Sherpa hadronisation module
# \ingroup python sherpa
#
# Base Sherpa configuration parameters to define the hadronisation module

from EventModifiers.sherpa_cfi import sherpa as _sherpa

sherpa = _sherpa.clone(
    seed = 1000,
    fragmentation = 'Ahadic',  # AHADIC++ cluster model; 'Pythia8' for the Lund string model (if built in Sherpa)
    hadronDecays = True,  # unstable hadrons decayed by HADRONS++
    runCard = '',  # optional Sherpa YAML run card for the fragmentation/decays tuning
    # any Sherpa setting (command-line syntax) steered here is applied before the first event is hadronised,
    # e.g. 'AHADIC: {M_UP_DOWN: 0.3}', 'HADRON_DECAYS: {Max_Proper_Lifetime: 100}', or 'AHADIC:M_UP_DOWN:0.3'
    preConfiguration = (),
    processConfiguration = (),
)
