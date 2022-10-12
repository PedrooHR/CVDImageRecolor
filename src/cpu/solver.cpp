/*
 * solver.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "solver.h"
#include <omp.h>

solver::solver(grid *A, dataset *B) {
  Grid = A;
  Dataset = B;

  // Size informations
  numEpochs = 7;
  zeroNode = -1;
  Emodule = 0.1;
  Rmodule = 50;
  CurrentEpoch = 0;
  GridDimension = 2;
  numNodes = A->graph_size;
  numPixels = Dataset->Datasize;

  // Final Solution, each Pixel has 3 values [L, a, b]
  Projection = (float *)malloc(numPixels * 3 * sizeof(float));

  // System Matrix
  Taxons = (int *)malloc(numPixels * sizeof(int));
  ResultVector = (float *)malloc(numNodes * 2 * sizeof(float));
  Solution = (float *)malloc(numNodes * 2 * sizeof(float));
  SysMatrix = (float *)malloc(numNodes * numNodes * sizeof(float));

  // Define the rates of Elastic Map for each epoch
  EpochM = (float *)malloc(numEpochs * sizeof(float));
  float epo = 0.05;
  for (int i = 0; i < numEpochs; i++) {
    EpochM[i] = epo;
    epo /= (3 - i * 0.04);
  }
}

void solver::resetSysMatrix() {
  for (int i = 0; i < numNodes; i++)
    for (int j = 0; j < numNodes; j++)
      SysMatrix[posN(i, j, numNodes)] = 0;
}

void solver::preConstructSysMatrix() {
  resetSysMatrix();

  Emodule = EpochM[CurrentEpoch] *
            pow(Grid->numEdges, (2 - GridDimension) / GridDimension) * 10;
  Rmodule = EpochM[CurrentEpoch] *
            pow(Grid->numRibs, (2 - GridDimension) / GridDimension) * 100;
}

void solver::constructSysMatrix() {
  int *TaxonSizes = (int *)malloc(numNodes * sizeof(int));

  // Calculate the centroids of each taxon
  for (int i = 0; i < numNodes; i++) {
    ResultVector[pos2(i, 0)] = 0;
    ResultVector[pos2(i, 1)] = 0;
    TaxonSizes[i] = 0;
  }
  for (int i = 0; i < numPixels; i++) {
    ResultVector[pos2(Taxons[i], 0)] += Dataset->Datapoints[a(i)];
    ResultVector[pos2(Taxons[i], 1)] += Dataset->Datapoints[b(i)];
    TaxonSizes[Taxons[i]] += 1;
  }
  for (int i = 0; i < numNodes; i++) {
    ResultVector[pos2(i, 0)] /= (1.0 * numPixels);
    ResultVector[pos2(i, 1)] /= (1.0 * numPixels);
  }

  // Influence of the ribs on the Elastic Map
  for (int i = 0; i < Grid->numRibs; i++) {
    int N0 = Grid->Ribs[i * 3 + 0];
    int N1 = Grid->Ribs[i * 3 + 1];
    int N2 = Grid->Ribs[i * 3 + 2];

    SysMatrix[posN(N0, N0, numNodes)] += Rmodule * 4.0;
    SysMatrix[posN(N0, N1, numNodes)] -= Rmodule * 2.0;
    SysMatrix[posN(N0, N2, numNodes)] -= Rmodule * 2.0;
    SysMatrix[posN(N1, N0, numNodes)] -= Rmodule * 2.0;
    SysMatrix[posN(N2, N0, numNodes)] -= Rmodule * 2.0;
    SysMatrix[posN(N1, N1, numNodes)] += Rmodule;
    SysMatrix[posN(N1, N2, numNodes)] += Rmodule;
    SysMatrix[posN(N2, N1, numNodes)] += Rmodule;
    SysMatrix[posN(N2, N2, numNodes)] += Rmodule;
  }

  // Influence of the Edges on the Elastic Map
  for (int i = 0; i < Grid->numEdges; i++) {
    int N0 = Grid->Edges[i * 2 + 0];
    int N1 = Grid->Edges[i * 2 + 1];

    SysMatrix[posN(N0, N0, numNodes)] += Emodule;
    SysMatrix[posN(N1, N1, numNodes)] += Emodule;
    SysMatrix[posN(N0, N1, numNodes)] -= Emodule;
    SysMatrix[posN(N1, N0, numNodes)] -= Emodule;
  }

  // Influence of the Taxons on the Elastic Map
  for (int i = 0; i < numNodes; i++) {
    SysMatrix[posN(i, i, numNodes)] += TaxonSizes[i] / (1.0 * numPixels);
  }
}

void solver::solveLS() {
  // escalonador
  for (int k = 0; k < numNodes; k++) {
    for (int j = k + 1; j < numNodes; j++) {
      float multi =
          SysMatrix[posN(j, k, numNodes)] / SysMatrix[posN(k, k, numNodes)];
      for (int i = k; i < numNodes; i++)
        SysMatrix[posN(j, i, numNodes)] -=
            SysMatrix[posN(k, i, numNodes)] * multi;

      ResultVector[pos2(j, 0)] -= ResultVector[pos2(k, 0)] * multi;
      ResultVector[pos2(j, 1)] -= ResultVector[pos2(k, 1)] * multi;
    }
  }
  // resolve
  for (int k = numNodes - 1; k >= 0; k--) {
    double sum_x = 0;
    double sum_y = 0;
    for (int j = k + 1; j < numNodes; j++) {
      sum_x += SysMatrix[posN(k, j, numNodes)] * Solution[pos2(j, 0)];
      sum_y += SysMatrix[posN(k, j, numNodes)] * Solution[pos2(j, 1)];
    }
    Solution[pos2(k, 0)] =
        (ResultVector[pos2(k, 0)] - sum_x) / SysMatrix[posN(k, k, numNodes)];
    Solution[pos2(k, 1)] =
        (ResultVector[pos2(k, 1)] - sum_y) / SysMatrix[posN(k, k, numNodes)];
  }

  for (int i = 0; i < numNodes; i++) {
    Grid->Location[pX(i)] = Solution[pos2(i, 0)];
    Grid->Location[pY(i)] = Solution[pos2(i, 1)];
  }
}

void solver::calcTaxons() {
#pragma omp parallel for
  for (int i = 0; i < numPixels; i++) {
    float MinDist = (float)MAX_VALUE;
    int MinNodeRef;
    // For each pixel, which ElMap graph node is closest one?
    for (int j = 0; j < numNodes; j++) {
      float x1 = Dataset->Datapoints[a(i)];
      float y1 = Dataset->Datapoints[b(i)];
      float x2 = Grid->Location[pX(j)];
      float y2 = Grid->Location[pY(j)];
      float localDist = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2));
      if (localDist < MinDist) {
        MinDist = localDist;
        MinNodeRef = j;
      }
    }
    Taxons[i] = MinNodeRef;
  }
}

void solver::drawRecolored(const char *filepath) {
  cimg_library::CImg<float> image(Dataset->width, Dataset->height, 1, 3, 0);
  for (unsigned int i = 0; i < numPixels; i++) {
    int x = i % Dataset->width;
    int y = i / Dataset->width;
    float *RGB =
        getRGBColor(Projection[L(i)], Projection[a(i)], Projection[b(i)]);
    for (int j = 0; j < 3; j++) {
      RGB[j] = (RGB[j] > 1.0) ? 1.0 : RGB[j];
      RGB[j] = (RGB[j] < 0.0) ? 0.0 : RGB[j];
    }
    image(x, y, 0) = RGB[0] * 255.0;
    image(x, y, 1) = RGB[1] * 255.0;
    image(x, y, 2) = RGB[2] * 255.0;
  }
  image.save_jpeg(filepath, 100);
}

void solver::centerWhite() {
  float min_distance = (float)MAX_VALUE;
  int zeroNode = numNodes / 2;
  for (int i = 0; i < numPixels; i++) {
    float x1 = Dataset->Datapoints[a(i)];
    float y1 = Dataset->Datapoints[b(i)];
    float x2 = 0;
    float y2 = 0;
    float localZeroDist = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2));
    if (localZeroDist < min_distance) {
      zeroNode = Taxons[i];
      min_distance = localZeroDist;
    }
  }

  printf("Origin Node: %d \n", zeroNode);

  int size_positive = numNodes - zeroNode;
  int size_negative = zeroNode;
  float x_step = fabs((size_positive > size_negative)
                          ? (Grid->A[0] / (1.0 * size_positive))
                          : (-Grid->C[0] / (1.0 * size_negative)));

  OriginalMap = (float *)malloc(numNodes * 2 * sizeof(float));

  float x_start = 0;
  for (int i = zeroNode; i < numNodes; i++) {
    OriginalMap[pos2(i, 0)] = x_start;
    OriginalMap[pos2(i, 1)] = Grid->miAB * x_start;
    x_start += x_step;
  }

  x_start = 0 - x_step;
  for (int i = zeroNode - 1; i >= 0; i--) {
    OriginalMap[pos2(i, 0)] = x_start;
    OriginalMap[pos2(i, 1)] = Grid->miBC * x_start;
    x_start -= x_step;
  }
}

void solver::projectPoints() {
#pragma omp parallel for
  for (int i = 0; i < numPixels; i++) {
    int point = i;
    int node = Taxons[i];

#if CENTERWHITE
    Projection[L(point)] = Dataset->Datapoints[L(point)];
    Projection[a(point)] = OriginalMap[pos2(node, 0)];
    Projection[b(point)] = OriginalMap[pos2(node, 1)];
#else
    Projection[L(point)] = Dataset->Datapoints[L(point)];
    Projection[a(point)] = Grid->CVDPosition[pX(node)];
    Projection[b(point)] = Grid->CVDPosition[pY(node)];
#endif
  }
}
