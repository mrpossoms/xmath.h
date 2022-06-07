#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::vec<5> a, r, expected;
    float s = RAND_F;

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;

        expected[i] = a[i] * s;
    }

    r = a * s;

    for (int i = 5; i--;)
    {
        assert(fabs(expected[i] - r[i]) < 1e-8);
    }
}
