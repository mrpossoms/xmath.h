#include ".test.h"

#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)

TEST
{
	xmath::mat<4, 4> M = {
		{ 4, 0, 0, 0, },
		{ 0, 0, 2, 0, },
		{ 0, 1, 2, 0, },
		{ 1, 0, 0, 1, }
	};

	xmath::mat<4, 4> exp_inv = {
		{ .25, 0, 0, 0, },
		{ 0, -1, 1, 0, },
		{ 0, 0.5, 0, 0, },
		{ -0.25, 0, 0, 1, }
	};

	M.invert();

    for (size_t i = 4; i--;)
    for (size_t j = 4; j--;)
    {
        assert(fabs(exp_inv[i][j] - M[i][j]) < 1e-8);
    }

    return 0;
}