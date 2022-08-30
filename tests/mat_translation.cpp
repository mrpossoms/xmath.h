#include ".test.h"

#include "xmath.h"
#include <iostream>


TEST
{
    for (unsigned i = 0; i < 1000; i++)
    {
        auto p = rand_unit<3>();
        auto M = xmath::mat<4, 4>::translation(p);

        assert(p.is_near(M * xmath::vec<3>{0, 0, 0}, 1e-13));           
    }

    return 0;
}
