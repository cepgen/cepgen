#!/bin/sh

#--- proper gcc/gfortran version
export CVMFS_DIR=/cvmfs/sft.cern.ch/lcg/releases
source ${CVMFS_DIR}/gcc/6.2.0/x86_64-centos7/setup.sh
#--- Delphes linking
source ${CVMFS_DIR}/ROOT/6.18.04-b1762/x86_64-centos7-gcc62-opt/bin/thisroot.sh
export DELPHES_DIR=${CVMFS_DIR}/delphes/3.4.2-7776a/x86_64-centos7-gcc62-opt
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:${DELPHES_DIR}/include
#--- Pythia linking
export PYTHIA6_DIR=${CVMFS_DIR}/MCGenerators/pythia6/429.2-c4089/x86_64-centos7-gcc62-opt
export PYTHIA8_DIR=${CVMFS_DIR}/MCGenerators/pythia8/243-ac0f1/x86_64-centos7-gcc62-opt
export PYTHIA8DATA=${PYTHIA8_DIR}/share/Pythia8/xmldoc
#--- extra utilities
export PYTHONHOME=${CVMFS_DIR}/Python/2.7.15-075d4/x86_64-centos7-gcc62-opt
export GSL_DIR=${CVMFS_DIR}/GSL/2.5-32fc5/x86_64-centos7-gcc62-opt
export HEPMC_DIR=${CVMFS_DIR}/HepMC/2.06.09-0a23a/x86_64-centos7-gcc62-opt
export LHAPDF_PATH=${CVMFS_DIR}/MCGenerators/lhapdf/6.2.2-8a3e6/x86_64-centos7-gcc62-opt
export TBB_DIR=${CVMFS_DIR}/tbb/2019_U1-5939b/x86_64-centos7-gcc62-opt
export DAVIX_DIR=${CVMFS_DIR}/Davix/0.7.1-f7fe6/x86_64-centos7-gcc62-opt
export VDT_DIR=${CVMFS_DIR}/vdt/0.4.2-84b8c/x86_64-centos7-gcc62-opt

echo "Environment prepared for LXPLUS"
