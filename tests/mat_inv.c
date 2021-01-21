
#include ".test.h"

#include "xmath.h"

TEST
{

    double A[2][2] = 
		{ {1, 2},
		  {4, 11},
		};

    double A_inverted[2][2];
    
    double A_inverted_expected[2][2] =
		{ {3.66666666666, -0.666666666666},
		  {-1.3333333333, 0.33333333333},
		};

    mat_inv(2, 2, A, A_inverted);
    
    __MAT_EQ(2, 2, A_inverted, A_inverted_expected);
   	
    return 0;
}
