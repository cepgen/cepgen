##
# \file
# \ingroup python
#
# Collection of parameters for the steering of a Vegas integration algorithm

from Config.Integration.plain_cfi import plain


class VegasIntegrationMode:
    """Type of integration to be performed"""
    Stratified = -1
    ImportanceOnly = 0
    Importance = 1


vegas = plain.clone('Vegas',
    numFunctionCalls = 50000,
    treat = True,  # smoothing of the integrand
    chiSqCut = 1.2,
    # VEGAS-specific parameters
    iterations = 10,
    alpha = 1.5,
    mode = VegasIntegrationMode.Importance,
    verbose = -1,
    loggingOutput = 'cerr',
)
