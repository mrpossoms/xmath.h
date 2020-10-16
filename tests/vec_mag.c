#include ".test.h"

#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)

TEST
{
    float a[5], b[5], expected = 0;

    for (int i = 5; i--;)
    {
        a[i] = RAND_F;
        b[i] = RAND_F;
        expected += a[i] * b[i];
    }

    expected = sqrt(expected);

    float actual = vec_mag(5, a, b);

    assert(fabs(expected - actual) < 1e-8);
}
