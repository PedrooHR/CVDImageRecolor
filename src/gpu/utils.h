#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>
#include <stdlib.h>
#include <vector>

// Execution definitions
#define USEOPENMP 1 // Only works if using CPU
#define USEGPU 1

// Mathematical definitions
#define euler 2.71828182846
#define SIGMA 1.0
#define MI 0.0

// Max/Min definitions
#define Min2(x, y) x < y ? x : y
#define Min3(x, y, z) Min2(x, Min2(y, z))
#define Max2(x, y) x > y ? x : y
#define Max3(x, y, z) Max2(x, Max2(y, z))

// Returns the index of Matrix with 2 columns
#define pos2(i, j) i * 2 + j
// Returns the index of Matrix with 3 columns
#define pos3(i, j) i * 3 + j
// Returns the index of Matrix with N columns
#define posN(i, j, N) i * N + j



// Operações entre espaço de cores
float *getRGBColor(float L, float a, float b);
float *getLabColor(unsigned int sR, unsigned int sG, unsigned int sB);

#endif
