#!/bin/bash

#source /afs/atlas.umich.edu/home/rkwang/public/sw/setup_TopDrawer.sh
#--- proper gcc/gfortran/cmake version
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9.3.0/x86_64-centos7-gcc9-opt/setup.sh

export LCG_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_101
export ARCH=x86_64-centos7-gcc8-opt

export CMAKE_DIR=${LCG_DIR}/CMake/3.20.0/${ARCH}
export PATH=${CMAKE_DIR}/bin:${PATH}
export Boost_DIR=${LCG_DIR}/Boost/1.77.0/${ARCH}
export GSL_DIR=${LCG_DIR}/GSL/2.7/${ARCH}
export OPENBLAS_DIR=${LCG_DIR}/blas/0.3.17.openblas/${ARCH}
#--- Delphes linking
source ${LCG_DIR}/ROOT/6.24.06/${ARCH}/bin/thisroot.sh
export DELPHES_DIR=${LCG_DIR}/delphes/3.5.0/${ARCH}
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:${DELPHES_DIR}/include
#--- Pythia linking
export PYTHIA6_DIR=${LCG_DIR}/MCGenerators/pythia6/429.2/${ARCH}
export PYTHIA8_DIR=${LCG_DIR}/MCGenerators/pythia8/245/${ARCH}
export PYTHIA8DATA=${PYTHIA8_DIR}/share/Pythia8/xmldoc
#--- Python environment
export PYTHONHOME=${LCG_DIR}/Python/3.9.6/${ARCH}
export PATH=${PYTHONHOME}/bin:${PATH}
#--- extra utilities
export APFEL_DIR=${LCG_DIR}/MCGenerators/apfel/3.0.4/${ARCH}
export HEPMC_DIR=${LCG_DIR}/HepMC/2.06.11/${ARCH}
export HEPMC3_DIR=${LCG_DIR}/hepmc3/3.2.4/${ARCH}
export LHAPDF_PATH=${LCG_DIR}/MCGenerators/lhapdf/6.3.0/${ARCH}
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current:${LHAPDF_DATA_PATH}
export PATH=${LHAPDF_PATH}/bin:${PATH}
export PHOTOSPP_DIR=${LCG_DIR}/MCGenerators/photos++/3.64/${ARCH}
export TAUOLAPP_DIR=${LCG_DIR}/MCGenerators/tauola++/1.1.8/${ARCH}
export TBB_DIR=${LCG_DIR}/tbb/2020_U2/${ARCH}
export VDT_DIR=${LCG_DIR}/vdt/0.4.3/${ARCH}
export LD_LIBRARY_PATH=${TBB_DIR}/lib:${VDT_DIR}/lib:${LD_LIBRARY_PATH}
export nlohmann_json_DIR=${LCG_DIR}/jsonmcpp/3.9.1/${ARCH}
export YODA_DIR=${LCG_DIR}/MCGenerators/yoda/1.9.0/${ARCH}

echo "Environment prepared for LXPLUS"
export CEPGEN_LXPLUS_ENV=1

export PATH=$PWD/install/bin:$PATH
export LD_LIBRARY_PATH=$PWD/install/lib64:$LD_LIBRARY_PATH
export CEPGEN_PATH=$PWD/install/share/CepGen
