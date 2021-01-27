#include ".test.h"

#include "xmath.h"

TEST
{

    double A[2][2] = 
		{ {1, 0},
		  {0, 1},
		};

    double B[2][2] = 
		{ {1, 0},
		  {0, 1},
		};

    double C[2][2];	
    double D[2][2];
    bzero(D, sizeof(double) * 4);	
	
    mat_sub(2, 2, C, A, B);
    __MAT_EQ(2, 2, C, D);
   	
    return 0;
}
