#include ".test.h"

#include "xmath.h"
#include <random>
#include <fstream>
#include <sstream>

#include "ascii.h"

using namespace ascii;
using namespace xmath;


constexpr auto dt = 0.1;

mat<2, 2> stm_cart = {
	{ 1, dt },
	{ 0, 1  }
};

mat<2, 2> process_noise_covar_cart = {
	{ 0.01, 0 },
	{ 0, 0.01 },
};

mat<2, 1> control_to_state_cart = { { 0 }, {dt} };

mat<1, 2> state_to_measurement_cart = {
	{ 1, 0 },
};

mat<1, 1> measurement_noise_covar_cart = { 
	{ 0.1 }
};

std::default_random_engine rng;


void random_walk_velocity()
{

	std::normal_distribution<double> process_noise(0,0.01);
	std::normal_distribution<double> measurement_noise(0.0,0.1);
	std::normal_distribution<double> wander_noise(-1,1);

	vec<2> state = {}; // consists of 1D position, and velocity

	filter::kalman<1, 2, 1> cart_filter(state); 
	std::unordered_map<std::string, std::vector<double>> traces;

	int i = 0;

	state[1] = 1;

	for (float t = 0; t < 10; t += dt)
	{
		state[1] += wander_noise(rng);
		state = stm_cart * state + vec<2>{process_noise(rng),process_noise(rng)};

		if (i % 10 == 0)
		{
			auto z = state[0] + measurement_noise(rng);


			cart_filter.measurement_update(
				stm_cart,
				state_to_measurement_cart,
				mat<1, 1>{{z}},
				measurement_noise_covar_cart
			);
		}

		cart_filter.time_update(stm_cart, control_to_state_cart, {}, process_noise_covar_cart);

		traces["true pos"].push_back(state[0]);
		traces["true vel"].push_back(state[1]);
		traces["estimated pos"].push_back(cart_filter.estimated.state[0][0]);
		traces["estimated vel"].push_back(cart_filter.estimated.state[1][0]);
		i += 1;
	}

	Asciichart chart(traces);
	std::cerr << __FUNCTION__ << std::endl;
	std::cerr << chart.offset(3).height(30).legendPadding(3).Plot();

	// ensure each state element is within 1-sigma of the true state
	std::cerr << cart_filter.estimated.covariance.to_string() << std::endl;
	const auto& cov = cart_filter.estimated.covariance;
	assert(fabs(cart_filter.estimated.state[0][0] - state[0]) < sqrt(cov[0][0]));
	assert(fabs(cart_filter.estimated.state[1][0] - state[1]) < sqrt(cov[1][1]));
}

void constant_velocity()
{
	std::normal_distribution<double> process_noise(0,0.01);
	std::normal_distribution<double> measurement_noise(0.0,0.1);

	vec<2> state = {}; // consists of 1D position, and velocity

	filter::kalman<1, 2, 1> cart_filter(state); 
	std::unordered_map<std::string, std::vector<double>> traces;

	int i = 0;

	state[1] = 1;

	for (float t = 0; t < 10; t += dt)
	{
		state = stm_cart * state + vec<2>{process_noise(rng),process_noise(rng)};

		if (i % 10 == 0)
		{
			auto z = state[0] + measurement_noise(rng);

			cart_filter.measurement_update(
				stm_cart,
				state_to_measurement_cart,
				mat<1, 1>{{z}},
				measurement_noise_covar_cart
			);
		}

		cart_filter.time_update(stm_cart, control_to_state_cart, {}, process_noise_covar_cart);

		traces["true pos"].push_back(state[0]);
		traces["true vel"].push_back(state[1]);
		traces["estimated pos"].push_back(cart_filter.estimated.state[0][0]);
		traces["estimated vel"].push_back(cart_filter.estimated.state[1][0]);
		i += 1;
	}

	Asciichart chart(traces);
	std::cerr << __FUNCTION__ << std::endl;
	std::cerr << chart.offset(3).height(30).legendPadding(10).Plot();

	// ensure each state element is within 1-sigma of the true state
	std::cerr << cart_filter.estimated.covariance.to_string() << std::endl;
	const auto& cov = cart_filter.estimated.covariance;
	assert(fabs(cart_filter.estimated.state[0][0] - state[0]) < sqrt(cov[0][0]));
	assert(fabs(cart_filter.estimated.state[1][0] - state[1]) < sqrt(cov[1][1]));
}

void sin_velocity()
{
	std::normal_distribution<double> process_noise(0,0.01);
	std::normal_distribution<double> measurement_noise(0.0,0.1);

	vec<2> state = {}; // consists of 1D position, and velocity

	filter::kalman<1, 2, 1> cart_filter(state); 
	std::unordered_map<std::string, std::vector<double>> traces;

	int i = 0;

	state[1] = 1;

	for (float t = 0; t < 10; t += dt)
	{
		state = stm_cart * state + vec<2>{process_noise(rng),process_noise(rng)};
		state[1] = sin(t);

		if (i % 1 == 0)
		{
			auto z = state[0] + measurement_noise(rng);

			cart_filter.measurement_update(
				stm_cart,
				state_to_measurement_cart,
				mat<1, 1>{{z}},
				measurement_noise_covar_cart
			);
		}

		cart_filter.time_update(stm_cart, control_to_state_cart, {}, process_noise_covar_cart);

		traces["true pos"].push_back(state[0]);
		traces["true vel"].push_back(state[1]);
		traces["estimated pos"].push_back(cart_filter.estimated.state[0][0]);
		traces["estimated vel"].push_back(cart_filter.estimated.state[1][0]);

		i += 1;
	}

	Asciichart chart(traces);
	std::cerr << __FUNCTION__ << std::endl;
	std::cerr << chart.offset(3).height(30).legendPadding(3).Plot();

	// ensure each state element is within 1-sigma of the true state
	std::cerr << cart_filter.estimated.covariance.to_string() << std::endl;
	const auto& cov = cart_filter.estimated.covariance;
	assert(fabs(cart_filter.estimated.state[0][0] - state[0]) < sqrt(cov[0][0]));
	assert(fabs(cart_filter.estimated.state[1][0] - state[1]) < sqrt(cov[1][1]));
}


TEST
{
	constant_velocity();
	random_walk_velocity();
	sin_velocity();
	return 0;
}
