/*
 * dataset.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "dataset.h"

dataset::dataset(const char *imgPath) {
  int local_id = 0;
  cimg_library::CImg<unsigned int> image(imgPath);
  width = image.width();
  height = image.height();
  int size = width * height + 1;

  Datapoints = (float *)malloc(size * 3 * sizeof(float));

  cimg_forXY(image, x, y) {
    // create Lab
    float *lab_color =
        getLabColor(image(x, y, 0), image(x, y, 1), image(x, y, 2));

    // Fill data points [0: L; 1: a; 2: b]
    Datapoints[L(local_id)] = lab_color[0];
    Datapoints[a(local_id)] = lab_color[1];
    Datapoints[b(local_id)] = lab_color[2];

    local_id++;
  };
  Datasize = size;
}
