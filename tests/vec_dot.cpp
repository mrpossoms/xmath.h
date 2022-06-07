#include ".test.h"

#include "xmath.h"

TEST
{
    using namespace xmath;

    { // check dot
        vec<5> a, b;
        float expected = 0;

        for (int i = 5; i--;)
        {
            a[i] = RAND_F;
            b[i] = RAND_F;
            expected += a[i] * b[i];
        }

        float actual = a.dot(b);

        assert(fabs(expected - actual) < 1e-8);   
    }

    { // check angle
        assert(fabs(vec<2>{1, 0}.angle_to(vec<2>{0, 1}) - (M_PI / 2)) < 1e-8);

        assert(fabs(vec<2>{1, 0}.angle_to(vec<2>{1, 0}) - 0) < 1e-8);

        assert(fabs(vec<2>{1, 0}.angle_to(vec<2>{-1, 0}) - M_PI) < 1e-8);
    }
}
