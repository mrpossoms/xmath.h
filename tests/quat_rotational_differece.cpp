#include ".test.h"
#include <cstddef>
#include <iostream>
#include "xmath.h"

using namespace xmath;

#define RAND_F (((random() % 2048) / 1024.f) - 1.f)

vec<3> rand_unit() 
{
    return vec<3>{ RAND_F, RAND_F, RAND_F}.unit();
}

TEST
{
    auto q0 = quat<>::from_axis_angle({1, 0, 0}, M_PI / 2);
    auto q1 = quat<>::from_axis_angle({1, 0, 0}, M_PI);

    assert(std::abs(q0.rotational_difference(q0)) < 0.0001);
    assert(std::abs(q0.rotational_difference(q1) - (M_PI / 2)) < 0.0001);

    // fuzz
    for (unsigned i = 0; i < 100; i--)
    {
        auto axis = rand_unit();
        auto a0 = RAND_F + 1;
        auto a1 = a0 + RAND_F + 1;

        q0 = quat<>::from_axis_angle(axis, a0);
        q1 = quat<>::from_axis_angle(axis, a1);
       
        assert(std::abs(q0.rotational_difference(q1) - (a1 - a0)) < 0.0001);
    }

    return 0;
}
