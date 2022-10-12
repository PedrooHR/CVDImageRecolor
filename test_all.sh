#!/bin/bash

g++ -c -pg src/solver.cpp
g++ -c -pg src/dataset.cpp
g++ -c -pg src/utils.cpp
g++ -c -pg src/grid.cpp
g++ -c -pg src/main.cpp
g++ *.o -pg -o main -lX11 -lpthread

rm -rf results
mkdir results

./main imgs/testeimg1.jpg
mv svimg.jpg results/testeimg1_a.jpg

./main imgs/testeimg2.jpg
mv svimg.jpg results/testeimg2_a.jpg

./main imgs/testeimg3.jpg
mv svimg.jpg results/testeimg3_a.jpg

./main imgs/testeimg4.jpg
mv svimg.jpg results/testeimg4_a.jpg

./main imgs/testeimg5.jpg
mv svimg.jpg results/testeimg5_a.jpg

./main imgs/testeimg6.jpg
mv svimg.jpg results/testeimg6_a.jpg

./main imgs/testeimg7.jpg
mv svimg.jpg results/testeimg7_a.jpg

./main imgs/testeimg8.jpg
mv svimg.jpg results/testeimg8_a.jpg

./main imgs/testeimg9.jpg
mv svimg.jpg results/testeimg9_a.jpg

./main imgs/testeimg10.jpg
mv svimg.jpg results/testeimg10_a.jpg

./main imgs/testeimg11.jpg
mv svimg.jpg results/testeimg11_a.jpg

./main imgs/testeimg12.jpg
mv svimg.jpg results/testeimg12_a.jpg

./main imgs/testeimg13.jpg
mv svimg.jpg results/testeimg13_a.jpg

./main testeimg14.jpg
mv svimg.jpg results/testeimg14_a.jpg

./main imgs/testeimg15.jpg
mv svimg.jpg results/testeimg15_a.jpg

./main imgs/testeimg16.jpg
mv svimg.jpg results/testeimg16_a.jpg

./main imgs/irina.jpg
mv svimg.jpg results/irina_a.jpg

exit
