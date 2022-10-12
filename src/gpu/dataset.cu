/*
 * dataset.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "dataset.h"

dataset::dataset(const char *imgPath) {
  cimg_library::CImg<unsigned int> image(imgPath);
  width = image.width();
  height = image.height();
  int size = width * height;
  
  cudaMallocManaged(&Datapoints, size * 3 * sizeof(float));

  for (int i = 0; i < size; i++) {
    int x = i % width;
    int y = i / width;
    // create Lab
    float *lab_color =
        getLabColor(image(x, y, 0), image(x, y, 1), image(x, y, 2));

    // Fill data points [0: L; 1: a; 2: b]
    Datapoints[L(i)] = lab_color[0];
    Datapoints[a(i)] = lab_color[1];
    Datapoints[b(i)] = lab_color[2];
  };
  Datasize = size;
}
