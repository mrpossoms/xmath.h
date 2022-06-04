#include ".test.h"
#include <cstddef>
#include <iostream>
#include "xmath.h"

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)


TEST
{
    using namespace xmath;

    { // Check simple, no rotation view
        auto q = quat<>::view({0, 0, 1}, {0, 1, 0});
        assert(q.is_near(quat<>{}, 0.00001f));
    }

    { // Check rotation in side plane
        auto q = quat<>::view({0, 1, 0}, {0, 0, -1});
        auto rotated = q.rotate(vec<3>{0, 0, 1});
        std::cerr << rotated.to_string() << std::endl;
        assert(rotated.is_near(vec<3>{0, 1, 0}, 0.00001f));
    }

    for (unsigned i = 100; i--;)
    {
        vec<3> target = { (RAND_F) * 10, (RAND_F) * 10, (RAND_F) * 10 };

        auto forward = target.unit();

        auto left = vec<3>::cross(forward, {0, -1, 0});
        auto q = quat<>::view(forward, vec<3>::cross(forward, left));

        auto rotated = q.rotate({0, 0, 1});

        auto delta = rotated - forward;
        auto mag = delta.magnitude();

        if (mag < 0.0001f)
        {
            std::cerr << "expected: " << forward.to_string() << " actual: " << rotated.to_string() << std::endl;
            std::cerr << mag << std::endl;
            assert(false);
        }
    }

    return 0;
}
