##
# \file
# \ingroup python
#
# Enumeration of the types of process kinematics (elastic, dissociative emissions)


class ProcessMode:
    """Types of processes supported"""
    ElectronProton      = 0
    ElasticElastic      = 1
    ElasticInelastic    = 2
    InelasticElastic    = 3
    InelasticInelastic  = 4
