#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::mat<2, 4> M = {
        { 1, 2, 3, 4, },
        { 5, 6, 7, 8, },
    };

    xmath::mat<2, 4> exp = {
        { 0, 0, 0, 0 },
        { 0, 0, 0, 0 },
    };

    auto O = M - M;

    for (size_t i = 2; i--;)
    for (size_t j = 4; j--;)
    {
        assert(fabs(exp[i][j] - O[i][j]) < 1e-8);
    }

    return 0;
}
