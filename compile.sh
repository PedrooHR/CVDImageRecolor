#!/bin/bash

g++ -c -pg src/solver.cpp
g++ -c -pg src/dataset.cpp 
g++ -c -pg src/utils.cpp
g++ -c -pg src/grid.cpp
g++ -c -pg src/main.cpp
g++ *.o -pg -o main -lX11 -lpthread -ljpeg 