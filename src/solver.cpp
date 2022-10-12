/*
 * solver.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "solver.h"

#define Min2(x, y) x < y ? x : y
#define Min3(x, y, z) Min2(x, Min2(y, z))
#define Max2(x, y) x > y ? x : y
#define Max3(x, y, z) Max2(x, Max2(y, z))

solver::solver(grid *A, dataset *B) {
  inverte = false;
  Grid = A;
  Dataset = B;
  numEpochs = 7;
  zeroNode = -1;
  Emodule = 0.1;
  Rmodule = 50;
  CurrentEpoch = 0;
  GridDimension = 2;
  numNodes = A->GRAPH_SIZE;
  SysMatrix.resize(numNodes);
  Taxons.resize(numNodes);
  ResultVector.resize(numNodes);
  Solution.resize(numNodes);
  NodesSolution.resize(Dataset->Datasize);
  for (int i = 0; i < numNodes; i++) {
    SysMatrix[i].resize(numNodes);
    Solution[i].resize(2);
    ResultVector[i].resize(2);
  }
  float epo = 0.05;
  for (int i = 0; i < numEpochs; i++) {
    EpochM.push_back(epo);
    epo /= (3 - i * 0.04);
  }
}

void solver::resetSysMatrix() {
  // Clear stage
  SysMatrix.clear();
  Taxons.clear();
  ResultVector.clear();
  Solution.clear();

  // Resize stage
  SysMatrix.resize(numNodes);
  Taxons.resize(numNodes);
  ResultVector.resize(numNodes);
  Solution.resize(numNodes);
  for (int i = 0; i < numNodes; i++) {
    SysMatrix[i].resize(numNodes);
    Solution[i].resize(2);
    ResultVector[i].resize(2);
  }
}

void solver::constructSysMatrix() {
  resetSysMatrix();

  Emodule = EpochM[CurrentEpoch] *
            pow(Grid->numEdges, (2 - GridDimension) / GridDimension) * 10;
  Rmodule = EpochM[CurrentEpoch] *
            pow(Grid->numRibs, (2 - GridDimension) / GridDimension) * 100;

  for (int i = 0; i < numNodes; i++)
    Taxons[i].clear();
  calcTaxons();

  // Funcional
  for (int i = 0; i < Grid->numRibs; i++) {
    int N0 = Grid->Ribs[i][0];
    int N1 = Grid->Ribs[i][1];
    int N2 = Grid->Ribs[i][2];

    SysMatrix[N0][N0] += Rmodule * 4.0;
    SysMatrix[N0][N1] -= Rmodule * 2.0;
    SysMatrix[N0][N2] -= Rmodule * 2.0;
    SysMatrix[N1][N0] -= Rmodule * 2.0;
    SysMatrix[N2][N0] -= Rmodule * 2.0;
    SysMatrix[N1][N1] += Rmodule;
    SysMatrix[N1][N2] += Rmodule;
    SysMatrix[N2][N1] += Rmodule;
    SysMatrix[N2][N2] += Rmodule;
  }

  for (int i = 0; i < Grid->numEdges; i++) {
    int N0 = Grid->Edges[i][0];
    int N1 = Grid->Edges[i][1];

    SysMatrix[N0][N0] += Emodule;
    SysMatrix[N1][N1] += Emodule;
    SysMatrix[N0][N1] -= Emodule;
    SysMatrix[N1][N0] -= Emodule;
  }

  for (int i = 0; i < numNodes; i++) {
    SysMatrix[i][i] += Taxons[i].size() / (1.0 * Dataset->Datasize);
  }

  for (int i = 0; i < numNodes; i++) {
    ResultVector[i][0] = 0;
    ResultVector[i][1] = 0;

    for (int j = 0; j < Taxons[i].size(); j++) {
      ResultVector[i][0] += Dataset->Datapoints[Taxons[i][j]].Lab[1];
      ResultVector[i][1] += Dataset->Datapoints[Taxons[i][j]].Lab[2];
    }
    ResultVector[i][0] /= (1.0 * Dataset->Datasize);
    ResultVector[i][1] /= (1.0 * Dataset->Datasize);
  }
}

void solver::solveLS() {
  // escalonador
  for (int k = 0; k < numNodes; k++) {
    for (int j = k + 1; j < numNodes; j++) {
      float multi = SysMatrix[j][k] / SysMatrix[k][k];
      for (int i = k; i < numNodes; i++)
        SysMatrix[j][i] -= (SysMatrix[k][i] * multi);

      ResultVector[j][0] -= ResultVector[k][0] * multi;
      ResultVector[j][1] -= ResultVector[k][1] * multi;
    }
  }
  // resolve
  for (int k = numNodes - 1; k >= 0; k--) {
    double somax = 0;
    double somay = 0;
    for (int j = k + 1; j < numNodes; j++) {
      somax += SysMatrix[k][j] * Solution[j][0];
      somay += SysMatrix[k][j] * Solution[j][1];
    }
    Solution[k][0] = (ResultVector[k][0] - somax) / SysMatrix[k][k];
    Solution[k][1] = (ResultVector[k][1] - somay) / SysMatrix[k][k];
  }

  for (int i = 0; i < numNodes; i++) {
    Grid->Location[i].x = Grid->Graph[i].Position[0] = Solution[i][0];
    Grid->Location[i].y = Grid->Graph[i].Position[1] = Solution[i][1];
  }
}

/* GPU KERNEL */
/****************************************************************************** 
__global__ void calcTaxonI(Pixel *Datapoints, Node *Graph, int size) {
  int i = threadIdx.x + (blockDim.x * blockIdx.x);

  if (i < size) {
    float MinDist = (float)MAX_VALUE;
    int MinNodeRef;
    for (int j = 0; j < T_GRID; j++) {

      float localDist =
          (Datapoints[i].a - Graph[j].a) * (Datapoints[i].a - Graph[j].a);
      localDist +=
          (Datapoints[i].b - Graph[j].b) * (Datapoints[i].b - Graph[j].b);

      if (localDist < MinDist) {
        MinDist = localDist;
        MinNodeRef = j;
      }
    }
    Datapoints[i].closest_node = MinNodeRef;
  }
}
******************************************************************************/

void solver::calcTaxons() {
  /* This was a GPU KERNEL - Can it be a FPGA kernel? */
  for (int i = 0; i < Dataset->Datasize; i++) {
    float MinDist = (float)MAX_VALUE;
    int MinNodeRef;
    // For each pixel, which ElMap graph node is closest one?
    for (int j = 0; j < numNodes; j++) {
      float x1 = Dataset->Location[i].x;
      float y1 = Dataset->Location[i].y;
      float x2 = Grid->Location[j].x;
      float y2 = Grid->Location[j].y; 
      float localDist = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2));
      if (localDist < MinDist) {
        MinDist = localDist;
        MinNodeRef = j;
      }
    }
    Taxons[MinNodeRef].push_back(i);
  }
}

void solver::drawRecolored(const char *filepath) {
  cimg_library::CImg<float> image(Dataset->width, Dataset->height, 1, 3, 0);
  for (int i = 0; i < NodesSolution.size(); i++) {
    int x = i % Dataset->width;
    int y = i / Dataset->width;
    std::vector<float> RGBcolor = getRGBColor(
        NodesSolution[i][2], NodesSolution[i][0], NodesSolution[i][1]);
    for (int j = 0; j < 3; j++) {
      RGBcolor[j] = (RGBcolor[j] > 1.0) ? 1.0 : RGBcolor[j];
      RGBcolor[j] = (RGBcolor[j] < 0.0) ? 0.0 : RGBcolor[j];
    }
    image(x, y, 0) = RGBcolor[0] * 255.0;
    image(x, y, 1) = RGBcolor[1] * 255.0;
    image(x, y, 2) = RGBcolor[2] * 255.0;
  }
  image.save_jpeg(filepath, 100);
}

void solver::centerWhite() {
  float min_distance = (float)MAX_VALUE;
  int zeroNode = numNodes / 2;
  for (int i = 0; i < numNodes; i++) {
    for (int j = 0; j < Taxons[i].size(); j++) {
      int point = Taxons[i][j];
      float x1 = Dataset->Datapoints[point].a;
      float y1 = Dataset->Datapoints[point].b;
      float x2 = 0;
      float y2 = 0; 
      float localZeroDist = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2));
      if (localZeroDist < min_distance) {
        zeroNode = i;
        min_distance = localZeroDist;
      }
    }
  }
  printf("Origin Node: %d \n", zeroNode);

  int size_positive = numNodes - zeroNode;
  int size_negative = zeroNode;
  float x_step = fabs((size_positive > size_negative)
                          ? (Grid->A[0] / (1.0 * size_positive))
                          : (-Grid->C[0] / (1.0 * size_negative)));

  OriginalMap.resize(numNodes);

  float x_start = 0;
  for (int i = zeroNode; i < numNodes; i++) {
    std::vector<float> t = {x_start, Grid->miAB * x_start};
    OriginalMap[i] = t;
    x_start += x_step;
  }

  x_start = 0 - x_step;
  for (int i = zeroNode - 1; i >= 0; i--) {
    std::vector<float> t = {x_start, Grid->miBC * x_start};
    OriginalMap[i] = t;
    x_start -= x_step;
  }
}

void solver::projectPoints() {
  NodesSolution.clear();
  NodesSolution.resize(Dataset->Datasize);
  for (int i = 0; i < Taxons.size(); i++) { // para cada nodo
    for (int j = 0; j < Taxons[i].size(); j++) {
      int current_point = Taxons[i][j];
      std::vector<float> algproj, original;

      algproj = {Grid->Graph[i].CVDposition[0], Grid->Graph[i].CVDposition[1],
                 Dataset->Datapoints[current_point].Lab[0]};

      if (OriginalMap.size() > 0) {
        algproj = {OriginalMap[i][0], OriginalMap[i][1],
                   Dataset->Datapoints[current_point].Lab[0]};
      }

      NodesSolution[current_point].push_back(algproj[0]);
      NodesSolution[current_point].push_back(algproj[1]);
      NodesSolution[current_point].push_back(algproj[2]);
    }
  }
}
