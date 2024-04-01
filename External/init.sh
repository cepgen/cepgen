#!/bin/sh

PYTHIA_VERSION=6428
PYTHIA_DIR=./
PYTHIA_FILE=pythia${PYTHIA_VERSION}.f
if [ ! -f ${PYTHIA_DIR}${PYTHIA_FILE} ]; then
  curl "https://pythia.org/download/pythia6/${PYTHIA_FILE}" > ${PYTHIA_DIR}${PYTHIA_FILE}
fi

HERWIG_VERSION=6521
HERWIG_FILE=herwig${HERWIG_VERSION}.f
HERWIG_INC_FILE=HERWIG65.INC
HERWIG_INC2_FILE=herwig${HERWIG_VERSION}.INC
if [ ! -f ${HERWIG_FILE} ]; then
  curl "https://www.hep.phy.cam.ac.uk/theory/webber/Herwig/${HERWIG_FILE}" > ${HERWIG_FILE}
fi
if [ ! -f ${HERWIG_INC_FILE} ]; then
  curl "https://www.hep.phy.cam.ac.uk/theory/webber/Herwig/${HERWIG_INC_FILE}" > ${HERWIG_INC_FILE}
fi
if [ ! -f ${HERWIG_INC2_FILE} ]; then
  curl "https://www.hep.phy.cam.ac.uk/theory/webber/Herwig/${HERWIG_INC2_FILE}" > ${HERWIG_INC2_FILE}
fi

if [[ $(hostname -f) =~ ^lxplus[0-9]+.cern.ch ]]; then
  echo ">>> on LXPLUS"
  curl "https://raw.githubusercontent.com/BristolTopGroup/AnalysisSoftware/master/FindROOT.cmake" > cmake/FindROOT.cmake
fi
