/*
 * solver.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>

#include "dataset.h"
#include "grid.h"
#include "utils.h"

// Application Definitions
#define T_GRID 100
#define CENTERWHITE 1
#define MAX_VALUE 1e10

class solver {
public:
  grid *Grid;
  dataset *Dataset;
  int numNodes;
  int numEpochs;
  float MSE = MAX_VALUE;
  int GridDimension;
  int CurrentEpoch;
  float Rmodule;
  float Emodule;
  int zeroNode;
	int numPixels;

  float *OriginalMap;

  float *EpochM;

	// Definitions used in the System Matrix
  int *Taxons;
  float *SysMatrix;
  float *ResultVector;
  float *Solution;

	// Projection Result
  float *Projection;

  void drawRecolored(const char *filepath);

  solver(grid *A, dataset *B);
  void resetSysMatrix();
  void preConstructSysMatrix();
  void constructSysMatrix();
  void calcTaxons();
  void solveLS();
  void projectPoints();
  void centerWhite();
};

#endif /* SOLVER_H_ */
