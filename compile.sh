#!/bin/bash

if [ "$1" == "gpu" ]; then
  nvcc -c src/$1/*.cu
  g++ *.o -L/usr/local/cuda/lib64 -o main_$1 -lcuda -lcudart -lX11 -lpthread -ljpeg 
elif [ "$1" == "cpu" ]; then
  g++ -c -fopenmp -pg src/$1/*.cpp
  g++ *.o -fopenmp -pg -o main_$1 -lX11 -lpthread -ljpeg 
fi

rm -rf *.o
rm -rf *.out
rm -rf *.txt
rm -rf *.jpg
