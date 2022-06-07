#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::vec<5> a, b, r, expected;

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;
        b[i] = RAND_F;
        expected[i] = a[i] + b[i] * 5;
    }

    VEC_ADD_SCL(5, r.v, a.v, b.v, 5);

    for (int i = 5; i--;)
    {
        assert(fabs(expected[i] - r[i]) < 1e-8);
    }
}
