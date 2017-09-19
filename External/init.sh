#!/bin/sh

PYTHIA_VERSION=6428

PYTHIA_DIR=./
PYTHIA_FILE=pythia${PYTHIA_VERSION}.f
if [ ! -f ${PYTHIA_DIR}${PYTHIA_FILE} ]; then
  curl "http://home.thep.lu.se/~torbjorn/pythia6/${PYTHIA_FILE}" > ${PYTHIA_DIR}${PYTHIA_FILE}
fi

if [[ $(hostname -f) =~ ^lxplus[0-9]+.cern.ch ]]; then
  echo ">>> on LXPLUS"
  curl "https://raw.githubusercontent.com/BristolTopGroup/AnalysisSoftware/master/FindROOT.cmake" > cmake/FindROOT.cmake
fi
