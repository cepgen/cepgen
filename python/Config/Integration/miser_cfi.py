##
# \file
# \ingroup python
#
# Collection of parameters for the steering of a MISER integration algorithm

from Config.Integration.plain_cfi import plain


miser = plain.clone('MISER',
    # MISER-specific parameters
    estimateFraction = 0.1,
    alpha = 2.,
    dither = 0.,
)
