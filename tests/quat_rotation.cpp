#include ".test.h"
#include <cstddef>
#include <iostream>
#include "xmath.h"


TEST
{
    using namespace xmath;

    auto q = quat<>::from_axis_angle({1, 0, 0}, M_PI / 2);
    auto rotated = q.rotate({0, 0, 1});

    assert(rotated.is_near({0, -1, 0}, 0.000001f));

    return 0;
}
