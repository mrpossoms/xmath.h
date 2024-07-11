#include ".test.h"

#include "xmath.h"

TEST
{
    xmath::vec<3> d = { 1, -1, 0 };
    xmath::vec<3> n = { 0, 1, 0 };
    assert(xmath::vec<3>::reflect(d, n).is_near({1, 1, 0}));

    return 0;
}
