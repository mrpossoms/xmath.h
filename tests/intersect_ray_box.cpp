#include ".test.h"

#include "xmath.h"

TEST
{
    { // trivial hit
        xmath::vec<3> ray_o = { 0, 0, -3 };
        xmath::vec<3> ray_d = { 0, 0, 1 };

        xmath::vec<3> box_o = { 0, 0, 0 };
        xmath::vec<3> box_sides[] = {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
        };

        assert(!isnan(xmath::intersect::ray_box(ray_o, ray_d, box_o, box_sides)));
    }

    { // trivial miss
        xmath::vec<3> ray_o = { 0, 0, -3 };
        xmath::vec<3> ray_d = { 0, 0, -1 };

        xmath::vec<3> box_o = { 0, 0, 0 };
        xmath::vec<3> box_sides[] = {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
        };

        assert(isnan(xmath::intersect::ray_box(ray_o, ray_d, box_o, box_sides)));
    }

    { // translational miss
        xmath::vec<3> ray_o[4] = {
            { 1.1, 0, -3 },
            {-1.1, 0, -3 },
            { 0,  1.1, -3 },
            { 0, -1.1, -3 },
        };
        xmath::vec<3> ray_d = { 0, 0, 1 };
        xmath::vec<3> box_o = { 0, 0, 0 };
        xmath::vec<3> box_sides[] = {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
        };

        for (unsigned i = 0; i < 4; i++)
        {
            assert(isnan(xmath::intersect::ray_box(ray_o[i], ray_d, box_o, box_sides)));
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
        xmath::vec<3> box_o = { 0, 0, 0 };
        xmath::vec<3> box_sides[] = {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
        };

        for (unsigned i = 0; i < 4; i++)
        {
            assert(!isnan(xmath::intersect::ray_box(ray_o[i], ray_d, box_o, box_sides)));
        }
    }

    { // glancing miss
        /*     1 x (0, 0, 0)
             |   1
 (-1, 0, -1) |________
             \   |
              \  |
               \ |
                \|
                 o (0, 0, -3)

        */
        xmath::vec<3> ray_o = { 0, 0, -3 };
        xmath::vec<3> ray_d = {-1.1, 0,  2 };

        xmath::vec<3> box_o = { 0, 0, 0 };
        xmath::vec<3> box_sides[] = {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
        };

        assert(isnan(xmath::intersect::ray_box(ray_o, ray_d, box_o, box_sides)));
    }

    return 0;
}
