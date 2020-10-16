#include ".test.h"

#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)

TEST
{
    float a[5], b[5], r[5], expected[5];

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;
        b[i] = RAND_F;
        expected[i] = a[i] - b[i];
    }

    vec_sub(5, r, a, b);

    for (int i = 5; i--;)
    {
        assert(fabs(expected[i] - r[i]) < 1e-8);
    }
}
