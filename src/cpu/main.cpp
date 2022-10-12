
#include "dataset.h"
#include "grid.h"
#include "solver.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>


std::string imgpath = "testeimg1.jpg";
grid *Grid;
solver *Solver;
dataset *Dataset;

void Init();
void Recolor();

// Programa Principal
int main(int argc, char **argv) {
  if (argc > 1)
    imgpath = std::string(argv[1]);

  Init();
  Recolor();
}

void Init(void) {
  Grid = new grid(PROTANOPE, T_GRID);
  printf("Grid constructed\n");
  Dataset = new dataset(imgpath.c_str());
  printf("Dataset constructed\n");
  Solver = new solver(Grid, Dataset);
  printf("Solver constructed\n");
}

void Recolor() {
  int max_epochs = Solver->numEpochs;
  while (Solver->CurrentEpoch < max_epochs) {
    printf("Doing Epoch %d\n", Solver->CurrentEpoch);
    Solver->preConstructSysMatrix();
    Solver->calcTaxons();
    Solver->constructSysMatrix();
    Solver->solveLS();
    Solver->CurrentEpoch++;
  }
#if CENTERWHITE
  Solver->centerWhite();
#endif
  printf("projection step\n");
  Solver->projectPoints();
  Solver->drawRecolored("svimg.jpg");
  printf("image saved as svimg.jpg\n");
}
