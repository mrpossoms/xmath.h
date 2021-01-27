#include ".test.h"

#include "xmath.h"

TEST
{
	auto M = xmath::mat<4, 4>::I();

	xmath::mat<4, 4> I = {
		{ 1, 0, 0, 0, },
		{ 0, 1, 0, 0, },
		{ 0, 0, 1, 0, },
		{ 0, 0, 0, 1, }
	};

    for (size_t i = 4; i--;)
    for (size_t j = 4; j--;)
    {
        assert(fabs(I[i][j] - M[i][j]) < 1e-8);
    }

    return 0;
}
