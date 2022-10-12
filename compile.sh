#!/bin/bash

g++ -c -fopenmp -pg src/cpu/solver.cpp
g++ -c -fopenmp -pg src/cpu/dataset.cpp 
g++ -c -fopenmp -pg src/cpu/utils.cpp
g++ -c -fopenmp -pg src/cpu/grid.cpp
g++ -c -fopenmp -pg src/cpu/main.cpp
g++ *.o -fopenmp -pg -o main -lX11 -lpthread -ljpeg 