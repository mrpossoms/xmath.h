#include ".test.h"

#include "xmath.h"


TEST
{
    using namespace xmath;

    { // check trivial projection
        auto v = vec<3>{1, 0, 0};
        auto n = vec<3>{0, 1, 0};
        assert(v.project_onto_plane(n).is_near(v, 1e-8));   
    }

    { // check zero projection
        auto v = vec<3>{0, 100, 0};
        auto n = vec<3>{0, 1, 0};
        assert(v.project_onto_plane(n).is_near(vec<3>{0, 0, 0}, 1e-8));   
    }

    { // check zero projection
        auto v = vec<3>{1, 100, 0};
        auto n = vec<3>{0, 1, 0};
        assert(v.project_onto_plane(n).is_near(vec<3>{1, 0, 0}, 1e-8));   
    }
}
