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
    const float n = 1;
    const float f = 10;
    auto P = xmath::mat<4, 4>::perspective(n, f, M_PI / 2.f, 1.f);
    std::cerr << "Projection mat" << std::endl << P.to_string() << std::endl;

    for (float x = -1; x < 0 ; x += 0.1f)
    {
        auto p_1 = P * xmath::mat<4,1>{{1}, {0}, {n}, {1}};

        for (float z = n + 1; z < f; z += 1.f)
        {
            auto p0 = xmath::mat<4, 1>{{x}, {0}, {z}, {1}};
            std::cerr << " X " << std::endl;
            std::cerr << p0.to_string() << std::endl;
            auto x = P * p0;
            std::cerr << " -> " << std::endl << x.to_string() << std::endl;
            // More distant points must be closer to the optical axis of the camera
            auto p = x / x[3][0];

            {

                // recover p0 from non-normalized projection
                assert((P.invert() * x).flatten().is_near(p0.flatten(), 0.0001f));

                if (!(std::abs(p[0][0]) < std::abs(p_1[0][0])))
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
                // b =  -fpn/fmn # P[2][2]
                // z_c = z * b - 1
                // w_c = z * a

                // recover z from z'
                // z' = ((z * b) - 1) / (z * a)
                // z' * (z * a) = (z * b) - 1
                // (z * z' * a) - (z * b) = -1
                // z * ((z' * a) - b) = -1
                // z = -1 / ((z' * a) - b)

                auto a = P[3][2];
                auto b = P[2][2];

                // given the derivation above, ensure that projection yeilds expected
                auto exp_z_prime = ((p0[2][0] * b) - 1) / (p0[2][0] * a);
                std::cerr << "exp_z_prime: " << exp_z_prime << " actual_z': " << p[2][0] << std::endl;
                assert(std::abs(exp_z_prime - p[2][0]) < 1e-6);

                auto _z = -1 / ((p[2][0] * a) - b);

                // undo the perspective divide and projection
                auto p_w = xmath::mat<4, 1>{
                    {p[0][0] * (_z * a)},
                    {p[1][0] * (_z * a)},
                    {_z},
                    {1}//{_z * a}
                };

                if (!p_w.flatten().is_near(p0.flatten(), 1e-6))
                {
                    std::cerr << "-" << std::endl;
                    std::cerr << "expected p_w: \n" << p0.to_string() << " actual p_w: \n" << p_w.to_string() << std::endl; 
                    assert(false);
                }
                else
                {
                    std::cerr << "+";
                }

                // auto p0_p = P.invert() * p_w;

                // if (!p0_p.flatten().is_near(x.flatten(), 1e-6))
                // {
                //     std::cerr << "-" << std::endl;
                //     std::cerr << std::endl << "expected z: " << p0[2][0] << " actual z: " << _z << std::endl;
                //     std::cerr << "expected p0: " << x.to_string() << " actual p0_p: " << p0_p.to_string() << std::endl; 
                //     assert(false);
                // }
                // else
                // {
                //     std::cerr << "+";
                // }

            }

            // std::cerr << p0.to_string() << " - " << x.to_string() << " - " << (x / x[3]).to_string() << std::endl;
        
        }
    }

    return 0;
}
