#include ".test.h"

#include "xmath.h"

TEST
{
	auto d = xmath::vec<4>{ 1, 2, 3, 4 };
	auto M = xmath::mat<4, 4>::diagonal(d);

	xmath::mat<4, 4> D = {
		{ 1, 0, 0, 0, },
		{ 0, 2, 0, 0, },
		{ 0, 0, 3, 0, },
		{ 0, 0, 0, 4, }
	};

    for (size_t i = 4; i--;)
    for (size_t j = 4; j--;)
    {
        assert(fabs(D[i][j] - M[i][j]) < 1e-8);
    }

    return 0;
}
