#!/bin/sh

if [ ! -d mstw2008grids ]; then
  if [ ! -f mstw2008grids.tar.gz ]; then
    wget "https://mstwpdf.hepforge.org/code/mstw2008grids.tar.gz" ;
  fi
  tar xfz mstw2008grids.tar.gz
fi

