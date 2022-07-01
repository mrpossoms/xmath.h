#include ".test.h"

#include "xmath.h"
#include <iostream>

TEST
{
    auto M = xmath::mat<4, 4>::look({1, 33, 7}, {0, 0, 1}, {0, 1, 0});

    auto p = M.transpose() * xmath::vec<3>{0, 0, 0};

    std::cerr << p.to_string() << "\n";
    assert(p.is_near({-1, 33, 7}));

    return 0;
}
