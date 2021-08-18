#!/bin/sh

#--- proper gcc/gfortran/cmake version
export CVMFS_DIR=/cvmfs/sft.cern.ch/lcg/releases
source ${CVMFS_DIR}/clang/9.0.0-a1c77/x86_64-centos7-gcc9-opt/setup.sh
source ${CVMFS_DIR}/gcc/9.3.0-467e1/x86_64-centos7/setup.sh
export CMAKE_DIR=${CVMFS_DIR}/CMake/3.17.3-75516/x86_64-centos7-gcc9-opt
export PATH=${CMAKE_DIR}/bin:${PATH}
#--- Delphes linking
source ${CVMFS_DIR}/ROOT/v6.24.00-f4a14/x86_64-centos7-gcc9-opt/bin/thisroot.sh
export DELPHES_DIR=${CVMFS_DIR}/delphes/3.4.3pre09-604d4/x86_64-centos7-gcc9-opt
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:${DELPHES_DIR}/include
#--- Pythia linking
export PYTHIA6_DIR=${CVMFS_DIR}/MCGenerators/pythia6/429.2-c4089/x86_64-centos7-gcc62-opt
export PYTHIA8_DIR=${CVMFS_DIR}/MCGenerators/pythia8/243-ac0f1/x86_64-centos7-gcc62-opt
export PYTHIA8DATA=${PYTHIA8_DIR}/share/Pythia8/xmldoc
#--- extra utilities
export PYTHONHOME=${CVMFS_DIR}/Python/3.9.5-82945/x86_64-centos7-gcc9-opt
export PATH=${PYTHONHOME}/bin:${PATH}
export GSL_DIR=${CVMFS_DIR}/GSL/2.5-32fc5/x86_64-centos7-gcc62-opt
export HEPMC_DIR=${CVMFS_DIR}/HepMC/2.06.09-0a23a/x86_64-centos7-gcc62-opt
export LHAPDF_PATH=${CVMFS_DIR}/MCGenerators/lhapdf/6.2.2-8a3e6/x86_64-centos7-gcc62-opt
export TBB_DIR=${CVMFS_DIR}/tbb/2020_U2-daa7e/x86_64-centos7-gcc9-opt
export VDT_DIR=${CVMFS_DIR}/vdt/0.4.3-992df/x86_64-centos7-gcc9-opt
export LD_LIBRARY_PATH=${TBB_DIR}/lib:${VDT_DIR}/lib:${LD_LIBRARY_PATH}
export JSONMCPP_DIR=${CVMFS_DIR}/jsonmcpp/3.9.1-72770/x86_64-centos7-gcc9-opt
export TAUOLAPP_DIR=${CVMFS_DIR}/MCGenerators/tauola++/1.1.6-3fa0d/x86_64-centos7-gcc62-opt

echo "Environment prepared for LXPLUS"
export CEPGEN_LXPLUS_ENV=1
