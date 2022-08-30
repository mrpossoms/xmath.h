#include ".test.h"

using namespace xmath;

TEST
{

    { // Test basic translation of the origin
        auto M = mat<4, 4>::look({1, 33, 7}, {0, 0, 1}, {0, 1, 0});

        auto p = M * vec<3>{1, 33, 7};

        std::cerr << p.to_string() << std::endl;
        assert(p.is_near({0, 0, 0}));
    }

    { // Both viewer and subject away from the origin, subject is in front of the viewer
        auto subject = vec<3>{11, 10, 10};
        auto M = mat<4, 4>::look_at({10, 10, 10}, subject, {0, 1, 0});

        auto p = M * subject;
        std::cerr << p.to_string() << std::endl;

        assert(p.dot({0, 0, 1}) >= 1.0);
    }

    for (unsigned i = 0; i < 100; i++)
    { // test orientation component
        auto f = rand_unit<3>();
        auto u = vec<3>{0, 1, 0};
        auto r = vec<3>::cross(f, u);
        u = vec<3>::cross(r, f).unit();

        auto O = mat<4, 4>{
            { r[0], u[0], f[0], 0 },
            { r[1], u[1], f[1], 0 },
            { r[2], u[2], f[2], 0 },
            { 0,    0,    0,    1 }
        }.transpose();

        auto p = O * f;
        std::cerr << p.to_string() << std::endl;
        auto actual = p.dot({0, 0, 1});
        assert(actual >= 0.99999);    
    }


    for (unsigned i = 0; i < 100; i++)
    { // same as above, but fuzz many random configurations
        auto view_pos = rand_unit<3>() * (rand() % 1024);
        auto view_dir = rand_unit<3>();

        auto subject_pos = view_pos + view_dir * (rand() % 100 + 1);

        auto M = mat<4, 4>::look_at(view_pos, subject_pos, {0, 1, 0});
        auto p = M * subject_pos;
        auto actual = p.dot({0, 0, 1});

        assert(actual >= 1.0);    
    }

    return 0;
}
