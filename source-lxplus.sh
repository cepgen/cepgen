#!/bin/sh

source /cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0/x86_64-centos7/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.16.00-4d56b/x86_64-centos7-gcc62-opt/bin/thisroot.sh
export PYTHONHOME=/cvmfs/sft.cern.ch/lcg/releases/Python/2.7.15-075d4/x86_64-centos7-gcc62-opt
export PYTHIA8_DIR=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia8/240p1-ecd34/x86_64-centos7-gcc62-opt
export PYTHIA8DATA=${PYTHIA8_DIR}/share/Pythia8/xmldoc
echo "Environment prepared for LXPLUS"
