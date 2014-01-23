#!/bin/sh

#gsl-randist 0 100000 gaussian 20 | gsl-histogram 0 100 200 > histogram.dat
#awk '{print $1, $3 ; print $2, $3}' histogram.dat | graph -T svg -C > test.svg

gsl-randist 0 100000 gaussian 20 | gsl-histogram 0 100 200 | awk '{print $1, $3 ; print $2, $3}' | graph -T svg -C > test.svg
