#include ".test.h"
#include <iostream>
#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)

TEST
{

	/* 	Square matrix case 	*/
	xmath::mat<4, 4> M = {
		{ 4, 0, 0, 0, },
		{ 0, 0, 2, 0, },
		{ 0, 1, 2, 0, },
		{ 1, 0, 0, 1, }
	};

	auto transpose_square = M.transpose();

	for (auto i = 0; i < 4; ++i)
	{
	for (auto j = 0; j < 4; ++j)
	{
		assert(transpose_square[i][j] == M[j][i]);
	}
	}

	/* 	General matrix case 	*/
	xmath::mat<4, 3> N = {
		{ 4, 0, 0},
		{ 0, 0, 2},
		{ 0, 1, 2},
		{ 1, 0, 0}
	};

        auto transpose_general = N.transpose();
	
	for (auto i = 0; i < 3; ++i)
	{
	for (auto j = 0; j < 4; ++j)
	{
		assert(transpose_general[i][j] == N[j][i]);
	}
	}

    return 0;
}
