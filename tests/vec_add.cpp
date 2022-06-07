#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::vec<5> a, b, r, expected;

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;
        b[i] = RAND_F;
        expected[i] = a[i] + b[i];
    }

    //vec_add(5, r, a, b);
    r = a + b;

    for (int i = 5; i--;)
    {
        assert(fabs(expected[i] - r[i]) < 1e-8);
    }
}
