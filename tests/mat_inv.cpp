#include ".test.h"

#include "xmath.h"


TEST
{
	xmath::mat<4, 4> M = {
		{ 4, 0, 0, 0, },
		{ 0, 0, 2, 0, },
		{ 0, 1, 2, 0, },
		{ 1, 0, 0, 1, }
	};

	xmath::mat<4, 4> I = {
		{ 1, 0, 0, 0, },
		{ 0, 1, 0, 0, },
		{ 0, 0, 1, 0, },
		{ 0, 0, 0, 1, }
	};

	auto inv = M.invert();

    // Multiply M by its inverse to get the identity
	M = M * inv;

    for (size_t i = 4; i--;)
    for (size_t j = 4; j--;)
    {
        assert(fabs(I[i][j] - M[i][j]) < 1e-8);
    }

    return 0;
}
