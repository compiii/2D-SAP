#!/bin/bash

mpic++ -std=c++1z -Wall -Wextra -O2 -pedantic -c src/*.cpp src/1D/*.cpp src/2D/*.cpp -fopenmp
mpic++ -o 2D_SAP.run *.o -fopenmp
rm *.o