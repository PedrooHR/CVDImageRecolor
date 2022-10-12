/*
 * grid.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef GRID_H_
#define GRID_H_
#include "utils.h"
#include <vector>

#define PROTANOPE 1
#define DEUTERANOPE 0
#define TRITANOPE 0

// Macros for accessing XY like positions
#define pX(idx) idx * 2 + 0
#define pY(idx) idx * 2 + 1

typedef struct Node {
  int id;
  float Weight;
  float *CVDposition;
  float *Position;
  float a;
  float b;
} Node;

class grid {
public:
  // CVD plane limits
  std::vector<float> A;
  std::vector<float> C;

	// Size information
  int graph_size;
  int numRibs;
  int numEdges;

	// Angles of CVD planes
  float miAB;
  float miBC;

  // Components of the Elastic Map
  int *Edges;
  int *Ribs;

  // Position on CVD Plane
  float *CVDPosition;

	// This only holds temporary final solution of EM points
  float *Location;

  grid(int type, int size);
};

#endif /* GRID_H_ */
