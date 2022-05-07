#include ".test.h"

#include "xmath.h"

TEST
{
    auto M = xmath::mat<4, 4>::orthographic(-1, 1, -5, 5, 2, -2);

    auto p0 = M * xmath::vec<4>{-5, 0, 0, 1};
    assert(p0.is_near({-1, 0, 0, 1}));

    auto p1 = M * xmath::vec<4>{5, 0, 0, 1};
    assert(p1.is_near({1, 0, 0, 1}));

    auto p2 = M * xmath::vec<4>{0, -2, 0, 1};
    assert(p2.is_near({0, -1, 0, 1}));

    auto p3 = M * xmath::vec<4>{0, 2, 0, 1};
    assert(p3.is_near({0, 1, 0, 1}));

    return 0;
}
