#include ".test.h"

#include "xmath.h"
#include <iostream>

static xmath::mat<4, 4> perspective_linear(float n, float f)
{
    auto d = f - n;

    return {
        { 1,         0,          0,         0 },
        { 0,         1,          0,         0 },
        { 0,         0,        1/d,         0 },
        { 0,         0,         -1,         0 }
    };

    // z' =  2 * (z) / f
    // 
}

TEST
{
    float n = 1;
    float f = 10;

    auto P = xmath::mat<4, 4>::perspective(n, f, M_PI / 2.f, 1.f);
    // auto P = perspective_linear(n, 20);
    // auto s0 = xmath::vec<4>{0, 0, 10, 1};
    // auto t = P * s0;
    auto p_1 = P * xmath::vec<4>{1, 0, n, 1};

    std::cerr << "Projection mat" << std::endl << P.to_string() << std::endl;

    for (float z = n + 1; z < f; z += 1.f)
    {
        auto p0 = xmath::vec<4>{1, 0, z, 1};
        auto x = P * p0;
        // More distant points must be closer to the optical axis of the camera
        auto p = x / x[3];

        {

            // recover p0 from non-normalized projection
            assert((P.invert() * x).is_near(p0, 0.0001f));

            if (!(std::abs(p[0]) < std::abs(p_1[0])))
            {
                std::cerr << "p[0] was not less than p_1[0]" << std::endl;
                std::cerr << "p_1: " << p_1.to_string() << " p: " << p.to_string() << std::endl;
                assert(false);
            }

            p_1 = p;
        }

        // convert to ndc and recover
        {
            // a = -2*ftn/fmn # P[3][2]
            // b =  -fpn/fmn # P [2][2]
            // z_c = z * b - 1
            // w_c = z * a

            // recover z from z'
            // z' = ((z * b) - 1) / (z * a)
            // z' * (z * a) = (z * b) - 1
            // (z * z' * a) - (z * b) = -1
            // z * (z' * a - b) = -1
            // z = -1 / ((z' * a) - b)
            auto _z = -1 / ((p[2] * P[3][2]) - P[2][2]);

            // undo the perspective divide
            auto p_w = xmath::vec<4>{
                p[0] * (_z * P[3][2]),
                p[1] * (_z * P[3][2]),
                _z,
                _z * P[3][2]
            };

            auto p0_p = P.invert() * p_w;

            if (std::abs(_z - p0[2]) > 1e-6)
            {
                std::cerr << "-" << std::endl;
                std::cerr << std::endl << "expected z: " << p0[2] << " actual z: " << _z << std::endl;
                assert(false);
            }
            else
            {
                std::cerr << "+";
            }

        }

        // std::cerr << p0.to_string() << " - " << x.to_string() << " - " << (x / x[3]).to_string() << std::endl;
    
    }

    return 0;
}
