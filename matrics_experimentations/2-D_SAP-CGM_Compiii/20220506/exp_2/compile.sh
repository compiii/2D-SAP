#!/bin/bash

module load openmpi/intel-opa/gcc/64/1.10.4-hfi

mpic++ -std=c++1z -Wall -Wextra -O2 -pedantic -c src/*.cpp src/1D/*.cpp src/2D/*.cpp
mpic++ -o 2D_SAP.run *.o
rm *.o