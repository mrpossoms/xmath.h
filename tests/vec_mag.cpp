#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::vec<5> a;
    float expected = 0;

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;
        expected += a[i] * a[i];
    }

    expected = sqrt(expected);
    float actual = a.magnitude();

    assert(fabs(expected - actual) < 1e-8);
}
