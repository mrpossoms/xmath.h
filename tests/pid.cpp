#include ".test.h"

#include "xmath.h"
#include <random>
#include <fstream>
#include <sstream>

#include "ascii.h"

using namespace ascii;
using namespace xmath;

constexpr auto dt = 0.1;

void P_converge()
{
    vec<1> x0 = {0};
    vec<1> x = x0;
    vec<1> setpoint = {10};

    xmath::controller::pid<1> pid({0.1});
    std::unordered_map<std::string, std::vector<double>> traces;

    for (unsigned i = 0; i < 100; i++)
    {
        x += pid(setpoint - x, dt);
        traces["X"].push_back(x[0]);
        traces["S"].push_back(setpoint[0]);
    }

    Asciichart chart(traces);
    std::cerr << __FUNCTION__ << std::endl;
    std::cerr << chart.offset(3).height(30).legend_padding(3).Plot();
    assert(pid.state.error[0] < 0.01);
}

void PD_converge()
{
    vec<1> x0 = {0};
    vec<1> x = x0;
    vec<1> p;
    vec<1> setpoint = {10};

    xmath::controller::pid<1> pid({0.1}, {}, {0.1});
    std::unordered_map<std::string, std::vector<double>> traces;

    for (unsigned i = 0; i < 100; i++)
    {
        p += pid(setpoint - x, dt);
        x += p * dt;
        traces["X"].push_back(x[0]);
        traces["S"].push_back(setpoint[0]);
    }

    Asciichart chart(traces);
    std::cerr << __FUNCTION__ << std::endl;
    std::cerr << chart.offset(3).height(30).legend_padding(3).Plot();
    assert(pid.state.error[0] < 0.01);
}

void PID_converge()
{
    vec<1> x0 = {0};
    vec<1> x = x0;
    vec<1> p;
    vec<1> setpoint = {10};

    xmath::controller::pid<1> pid({0.2}, {0.05}, {0.9});
    std::unordered_map<std::string, std::vector<double>> traces;

    for (unsigned i = 0; i < 100; i++)
    {
        auto u = pid(setpoint - x, dt);

        p -= {10 * dt};
        p += u;
        x += p * dt;

        if (x[0] < 0)
        {
            x[0] = 0;
            p[0] = 0;
        }

        traces["X"].push_back(x[0]);
        traces["S"].push_back(setpoint[0]);
    }

    Asciichart chart(traces);
    std::cerr << __FUNCTION__ << std::endl;
    std::cerr << chart.offset(3).height(30).legend_padding(3).Plot();
    assert(pid.state.error[0] < 0.1);
}

TEST
{
    P_converge();
    PD_converge();
    PID_converge();

    return 0;
}
