/*
 * grid.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "grid.h"

#include <iostream>

grid::grid(int type, int size) {
  graph_size = size;
  numEdges = graph_size - 1;
  numRibs = graph_size - 2;

  // Each pixel has 2
  CVDPosition = (float *)malloc(graph_size * 2 * sizeof(float));
  cudaMallocManaged(&Location, graph_size * 2 * sizeof(float));

  // do graphnodes
  if (type == PROTANOPE) {
    // PROTANOPE LIMITS
    A = {8.648425, -73.086372, 56.664734};
    C = {-14.907598, 86.293831, 89.536812};
  }

  miAB = A[1] / A[0];
  miBC = C[1] / C[0];

  float x_step = fabs((A[0] - C[0]) / (graph_size * 1.0));
  float x_start = C[0];

  for (int i = 0; i < graph_size; i++) {
    Location[pX(i)] = (x_start <= 0) ? miBC * x_start : miAB * x_start;
    Location[pY(i)] = -x_start;
    CVDPosition[pX(i)] = x_start;
    CVDPosition[pY(i)] = (x_start <= 0) ? miBC * x_start : miAB * x_start;
    x_start += x_step;
  }

  // Create Elastic Map Edges - Each edge is composed by two parts
  Edges = (int *)malloc(numEdges * 2 * sizeof(int));
  for (int i = 0; i < numEdges; i++) {
    Edges[i * 2 + 0] = i;
    Edges[i * 2 + 1] = i + 1;
  }

  // Create Elastic Map Ribs - Each ribs is composed by three parts
  Ribs = (int *)malloc(numRibs * 3 * sizeof(int));
  for (int i = 0; i < numRibs; i++) {
    Ribs[i * 3 + 0] = i + 1;
    Ribs[i * 3 + 1] = i;
    Ribs[i * 3 + 2] = i + 2;
  }
}