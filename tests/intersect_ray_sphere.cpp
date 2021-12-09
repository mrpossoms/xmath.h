#include ".test.h"

#include "xmath.h"

TEST
{
    { // trivial hit
        xmath::vec<3> ray_o = { 0, 0, -3 };
        xmath::vec<3> ray_d = { 0, 0, 1 };

        xmath::vec<3> sphere_o = { 0, 0, 0 };

        assert(!isnan(xmath::intersect::ray_sphere(ray_o, ray_d, sphere_o, 1)));
    }

    { // trivial miss
        xmath::vec<3> ray_o = { 0, 0, -3 };
        xmath::vec<3> ray_d = { 0, 0, -1 };

        xmath::vec<3> sphere_o = { 0, 0, 0 };

        assert(isnan(xmath::intersect::ray_sphere(ray_o, ray_d, sphere_o, 1)));
    }

    { // translational miss
        xmath::vec<3> ray_o[4] = {
            { 1.1, 0, -3 },
            {-1.1, 0, -3 },
            { 0,  1.1, -3 },
            { 0, -1.1, -3 },
        };
        xmath::vec<3> ray_d = { 0, 0, 1 };
        xmath::vec<3> sphere_o = { 0, 0, 0 };

        for (unsigned i = 0; i < 4; i++)
        {
            assert(isnan(xmath::intersect::ray_sphere(ray_o[i], ray_d, sphere_o, 1)));
        }
    }


    { // translational miss
        xmath::vec<3> ray_o[4] = {
            { 1.1, 0, -3 },
            {-1.1, 0, -3 },
            { 0,  1.1, -3 },
            { 0, -1.1, -3 },
        };
        xmath::vec<3> ray_d = { 0, 0, 1 };
        xmath::vec<3> sphere_o = { 0, 0, 0 };

        for (unsigned i = 0; i < 4; i++)
        {
            assert(isnan(xmath::intersect::ray_sphere(ray_o[i], ray_d, sphere_o, 1)));
        }
    }

    { // translational hit
        xmath::vec<3> ray_o[4] = {
            { 0.999, 0, -3 },
            {-0.999, 0, -3 },
            { 0,  0.999, -3 },
            { 0, -0.999, -3 },
        };
        xmath::vec<3> ray_d = { 0, 0, 1 };
        xmath::vec<3> sphere_o = { 0, 0, 0 };

        for (unsigned i = 0; i < 4; i++)
        {
            assert(!isnan(xmath::intersect::ray_sphere(ray_o[i], ray_d, sphere_o, 1)));
        }
    }

    return 0;
}
