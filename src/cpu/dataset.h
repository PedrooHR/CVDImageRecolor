/*
 * dataset.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef DATASET_H_
#define DATASET_H_

#define cimg_display 0
#define cimg_use_jpeg

#include <CImg.h>
#include "utils.h"
#include <vector>

// Macros for accessing Lab values
#define L(idx) idx * 3 + 0
#define a(idx) idx * 3 + 1
#define b(idx) idx * 3 + 2

class dataset {
public:
  int width;
  int height;
  int Datasize;

  // This represents a N by 3 array. Each N-th element is [L, a, b]. Each
  // element is a Pixel in Lab
  float *Datapoints;

  dataset(const char *imgPath);
};

#endif /* DATASET_H_ */
