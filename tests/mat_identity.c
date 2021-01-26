
#include ".test.h"

#include "xmath.h"

TEST
{

    double A[2][2] = 
		{ {1, 2},
		  {4, 11},
		};

    double A_Identity_Expected[2][2] =
		{ {1, 0},
		  {0, 1},
		};
   
    mat_identity(2, 2, A);  
    __MAT_EQ(2, 2, A, A_Identity_Expected);
   	
    return 0;
}
