#!/bin/sh

PYTHIA_VERSION=6428

PYTHIA_DIR=./
PYTHIA_FILE=pythia${PYTHIA_VERSION}.f
if [ ! -f ${PYTHIA_DIR}${PYTHIA_FILE} ]; then
  curl "http://home.thep.lu.se/~torbjorn/pythia6/${PYTHIA_FILE}" > ${PYTHIA_DIR}${PYTHIA_FILE}
fi
