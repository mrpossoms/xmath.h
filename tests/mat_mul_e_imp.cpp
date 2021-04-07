#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::mat<2, 4> M = {
        { 1, 2, 3, 4, },
        { 5, 6, 7, 8, },
    };


    xmath::mat<2, 4> exp = {
        { 10, 20, 30, 40, },
        { 50, 60, 70, 80, },
    };

    M *= 10.0;

    for (size_t i = 2; i--;)
    for (size_t j = 3; j--;)
    {
        assert(fabs(exp[i][j] - M[i][j]) < 1e-8);
    }

    return 0;
}
