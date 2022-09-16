#include ".test.h"

#include "xmath.h"
#include <random>
#include <fstream>

#include "ascii.h"

using namespace ascii;
using namespace xmath;


void cart_constant_velocity()
{
	
	mat<2, 1> state = {}; // consists of 1D position, and velocity

	filter::kalman<1, 2, 1> cart_filter(state); 
	constexpr auto dt = 0.1;

	mat<2, 2> stm = {
		{ 1, dt },
		{ 0, 1  }
	};

	mat<2, 2> process_noise_covar = {
		{ 0.01, 0 },
		{ 0, 0.01 },
	};

	mat<2, 1> control_to_state = { { 0 }, {dt} };

	mat<1, 2> state_to_measurement = {
		{ 1, 0 },
	};

	mat<1, 1> measurement_noise_covar = { 
		{ 0.1 }
	};

	int i = 0;

	state[1][0] = 1;
	std::default_random_engine rng;
	std::normal_distribution<double> process_noise(0,0.01);
	std::normal_distribution<double> measurement_noise(0.0,0.1);

	std::ofstream file("kalman-cart-constant.csv");
	std::vector<double> filtered_pos_data, truth_pos_data;
	for (float t = 0; t < 10; t += dt)
	{
		state = stm * state + mat<2,1>{{process_noise(rng)},{process_noise(rng)}};

		if (i % 10 == 0)
		{
			auto z = state[0] + measurement_noise(rng);

			cart_filter.measurement_update(
				stm,
				state_to_measurement,
				mat<1, 1>{{z}},
				measurement_noise_covar
			);
		}

		cart_filter.time_update(stm, control_to_state, {}, process_noise_covar);

		filtered_pos_data.push_back(cart_filter.estimated.state[0][0]);
		truth_pos_data.push_back(state[0][0]);

		file << t << "," << state[0][0] << "," << cart_filter.estimated.state[0][0] << ","
		                 << state[1][0] << "," << cart_filter.estimated.state[1][0] << "\n";

		i += 1;
	}

	file.close();

	Asciichart chart({
		{"truth pos", truth_pos_data}, 
		{"estimated pos", filtered_pos_data}
	});

	std::cerr << chart.offset(3).legendPadding(3).Plot();

	// ensure each state element is within 1-sigma of the true state
	std::cerr << cart_filter.estimated.covariance.to_string() << std::endl;
	const auto& cov = cart_filter.estimated.covariance;
	assert(fabs(cart_filter.estimated.state[0][0] - state[0][0]) < sqrt(cov[0][0]));
	assert(fabs(cart_filter.estimated.state[1][0] - state[1][0]) < sqrt(cov[1][1]));

}

void cart_random_walk_velocity()
{
	
	mat<2, 1> state = {}; // consists of 1D position, and velocity

	filter::kalman<1, 2, 1> cart_filter(state); 
	constexpr auto dt = 0.1;

	mat<2, 2> stm = {
		{ 1, dt },
		{ 0, 1  }
	};

	mat<2, 2> process_noise_covar = {
		{ 0.01, 0 },
		{ 0, 0.01 },
	};

	mat<2, 1> control_to_state = { { 0 }, {dt} };

	mat<1, 2> state_to_measurement = {
		{ 1, 0 },
	};

	mat<1, 1> measurement_noise_covar = { 
		{ 0.1 }
	};

	int i = 0;

	state[1][0] = 1;
	std::default_random_engine rng;
	std::normal_distribution<double> process_noise(0,0.01);
	std::normal_distribution<double> measurement_noise(0.0,0.1);


	std::unordered_map<std::string, std::vector<double>> traces;
	std::ofstream file("kalman.csv");
	for (float t = 0; t < 10; t += dt)
	{
		state = stm * state + mat<2,1>{{process_noise(rng)},{process_noise(rng)}};
		state[1][0] += dt * process_noise(rng);

		if (i % 10 == 0)
		{
			auto z = state[0] + measurement_noise(rng);

			cart_filter.measurement_update(
				stm,
				state_to_measurement,
				mat<1, 1>{{z}},
				measurement_noise_covar
			);
		}

		cart_filter.time_update(stm, control_to_state, {}, process_noise_covar);

		file << t << "," << state[0][0] << "," << cart_filter.estimated.state[0][0] << ","
		                 << state[1][0] << "," << cart_filter.estimated.state[1][0] << "\n";

		traces["true pos"].push_back(state[0][0]);
		traces["true vel"].push_back(state[1][0]);
		traces["estimated pos"].push_back(cart_filter.estimated.state[0][0]);
		traces["estimated vel"].push_back(cart_filter.estimated.state[1][0]);
		i += 1;
	}

	file.close();

	Asciichart chart(traces);
	std::cerr << chart.offset(3).height(30).legendPadding(3).Plot();

	// ensure each state element is within 1-sigma of the true state
	std::cerr << cart_filter.estimated.covariance.to_string() << std::endl;
	const auto& cov = cart_filter.estimated.covariance;
	assert(fabs(cart_filter.estimated.state[0][0] - state[0][0]) < sqrt(cov[0][0]));
	assert(fabs(cart_filter.estimated.state[1][0] - state[1][0]) < sqrt(cov[1][1]));

}

TEST
{
	cart_constant_velocity();
	cart_random_walk_velocity();
	return 0;
}
