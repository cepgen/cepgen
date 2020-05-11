#!/bin/sh

#--- proper gcc/gfortran version
source /cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0/x86_64-centos7/setup.sh
#--- Delphes linking
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.18.04-b1762/x86_64-centos7-gcc62-opt/bin/thisroot.sh
export DELPHES_DIR=/cvmfs/sft.cern.ch/lcg/releases/delphes/3.4.2-7776a/x86_64-centos7-gcc62-opt
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:${DELPHES_DIR}/include
#--- Pythia 8 linking
export PYTHIA8_DIR=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia8/243-ac0f1/x86_64-centos7-gcc62-opt
export PYTHIA8DATA=${PYTHIA8_DIR}/share/Pythia8/xmldoc
#--- extra utilities
export PYTHONHOME=/cvmfs/sft.cern.ch/lcg/releases/Python/2.7.15-075d4/x86_64-centos7-gcc62-opt

echo "Environment prepared for LXPLUS"
