#include ".test.h"
#include <cstddef>
#include <iostream>
#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)


TEST
{
    using namespace xmath;

    for (unsigned i = 100; i--;)
    {
        vec<3> target = { (RAND_F) * 10, (RAND_F) * 10, (RAND_F) * 10 };

        auto forward = target.unit();

        auto left = vec<3>::cross(forward, {0, -1, 0});
        auto q = quat<>::view(forward, vec<3>::cross(forward, left));

        auto rotated = q.rotate({0, 0, 1});

        auto delta = rotated - forward;
        auto mag = delta.magnitude();
        assert(mag < 0.0001f);

    }

    return 0;
}
