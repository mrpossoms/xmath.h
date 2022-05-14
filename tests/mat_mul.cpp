#include ".test.h"

#include "xmath.h"

TEST
{
    { // test matrix multiplication
        xmath::mat<2, 4> M = {
            { 1, 2, 3, 4, },
            { 5, 6, 7, 8, },
        };

        xmath::mat<4, 3> N = {
            { 1, 2, 3, },
            { 4, 5, 6, },
            { 7, 8, 9, },
            { 10,11,12,},
        };

        xmath::mat<2, 3> exp = {
            { 70, 80, 90, },
            { 158, 184, 210, },
        };

        auto O = M * N;

        for (size_t i = 2; i--;)
        for (size_t j = 3; j--;)
        {
            assert(fabs(exp[i][j] - O[i][j]) < 1e-8);
        }
    }

    { // test specialized matrix mult with vec<3>
        xmath::mat<4, 4> M = {
            { 1, 0, 0, 1 },
            { 0, 2, 0, 0 },
            { 0, 0, 1, 0 },
            { 0, 0, 0, 1 }
        };

        auto res = M * xmath::vec<3>{ 0, 1, 0 };

        xmath::vec<3> exp = {1, 2, 0};

        for (size_t j = 3; j--;)
        {
            assert(fabs(exp[j] - res[j]) < 1e-8);
        }   
    }

    return 0;
}
