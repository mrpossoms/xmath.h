#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <xmath.h>

#define TEST int main (int argc, const char* argv[])
#define TOLERANCE (1e-8)

#define RAND_F (((rand() % 2048) / 1024.f) - 1.f)

#define __MAT_EQ(rows, cols, left_matrix, right_matrix)\
{\
    for (size_t __i = (rows); __i--;)\
    for (size_t __j = (cols); __j--;)\
    {\
        assert(fabs(left_matrix[__i][__j] - right_matrix[__i][__j]) < TOLERANCE);\
    }\
}\

#ifdef __cplusplus
#include <iostream>

template<size_t N>
xmath::vec<N> rand_unit()
{
    xmath::vec<N> v;

    for (unsigned i = 0; i < N; i++)
    {
        v[i] = RAND_F;
    }

    return v.unit();
}

#endif
// -Wall doesn't like ending a file with a backlash
