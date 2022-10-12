/*
 * solver.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>
#include <iostream>

#include "grid.h"
#include "dataset.h"
#include "utils.h"

#define T_GRID 100
#define GPU_ENABLE 0

#define MAX_VALUE 1e10

class solver {

public:
	grid *Grid;
	dataset *Dataset;
	int numNodes;
	bool inverte;
	int numEpochs;
	float MSE = MAX_VALUE;
	int GridDimension;
	int CurrentEpoch;
	float Rmodule;
	float Emodule;
	int zeroNode;
	//Projeção final
	std::vector <std::vector<float> > OriginalMap;
	
	std::vector <std::vector <int> > Taxons;
	std::vector <float> EpochM;
	std::vector <std::vector <float> > SysMatrix;
	std::vector <std::vector <float> > ResultVector;
	std::vector <std::vector <float> > Solution;
	std::vector <std::vector <float> > NodesSolution;
	void drawRecolored(const char * filepath);

	solver(grid *A, dataset *B);
	void resetSysMatrix();
	void constructSysMatrix();
	void calcTaxons();
	void solveLS();
	void projectPoints();
	void centerWhite();
};

#endif /* SOLVER_H_ */
