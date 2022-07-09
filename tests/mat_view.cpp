#include ".test.h"

TEST
{

    // { // Test basic translation of the origin
    //     auto M = xmath::mat<4, 4>::look({1, 33, 7}, {0, 0, 1}, {0, 1, 0});

    //     auto p = M * xmath::vec<3>{0, 0, 0};

    //     std::cerr << p.to_string() << std::endl;
    //     assert(p.is_near({-1, -33, -7}));
    // }

    { // Both viewer and subject away from the origin, subject is in front of the viewer
        auto subject = xmath::vec<3>{11, 10, 10};
        auto M = xmath::mat<4, 4>::look_at({10, 10, 10}, subject, {0, 1, 0});

        auto p = M * subject;
        std::cerr << p.to_string() << std::endl;

        assert(p.dot({0, 0, 1}) >= 1.0);
    }

    for (unsigned i = 0; i < 100; i++)
    { // test orientation component
        auto f = rand_unit<3>();
        auto u = xmath::vec<3>{0, 1, 0};
        auto r = xmath::vec<3>::cross(f, u);
        u = xmath::vec<3>::cross(r, f).unit();

        auto O = xmath::mat<4, 4>{
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

        auto M = xmath::mat<4, 4>::look_at(view_pos, subject_pos, {0, 1, 0});
        auto p = M * subject_pos;

        auto actual = p.dot({0, 0, 1});
        std::cerr << p.to_string() << ": " << actual << std::endl;
        assert(actual >= 1.0);    
    }

    return 0;
}
