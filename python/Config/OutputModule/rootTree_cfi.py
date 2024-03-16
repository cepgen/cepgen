##
# \file
# \ingroup python root
#
# Collection of parameters to steer the ROOT TTree output module


from Config.containers_cff import Module


rootTree = Module('root_tree')

rootTreeCompressed = rootTree.clone(compressed = True)
