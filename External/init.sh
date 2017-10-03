#!/bin/sh

PYTHIA_VERSION=6428

PYTHIA_DIR=./
PYTHIA_FILE=pythia${PYTHIA_VERSION}.f
if [ ! -f ${PYTHIA_DIR}${PYTHIA_FILE} ]; then
  curl "http://home.thep.lu.se/~torbjorn/pythia6/${PYTHIA_FILE}" > ${PYTHIA_DIR}${PYTHIA_FILE}
fi

if [ ! -d mstw2008grids ]; then
  if [ ! -f mstw2008grids.tar.gz ]; then
    wget "https://mstwpdf.hepforge.org/code/mstw2008grids.tar.gz" ;
  fi
  tar xfz mstw2008grids.tar.gz
fi
