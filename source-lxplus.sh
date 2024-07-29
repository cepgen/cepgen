#!/bin/bash

#source /afs/atlas.umich.edu/home/rkwang/public/sw/setup_TopDrawer.sh
#--- proper gcc/gfortran/cmake version
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/13/x86_64-el9-gcc13-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

echo "Environment prepared for LXPLUS"
export CEPGEN_LXPLUS_ENV=1

export PATH=$PWD/install/bin:$PATH
export LD_LIBRARY_PATH=$PWD/install/lib64:$LD_LIBRARY_PATH
export CEPGEN_PATH=$PWD/install/share/CepGen
