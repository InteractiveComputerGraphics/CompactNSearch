
#include <CompactNSearch>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

#include <omp.h>

using namespace CompactNSearch;

std::vector<std::array<Real, 3>> positions;

std::size_t const N = 150;
Real const r_omega = 0.15;
Real const r_omega2 = r_omega * r_omega;
Real const radius = 2.0 * (2.0 * r_omega / static_cast<Real>(N-1));

std::size_t const N_enright_steps = 50;

Real
compute_average_number_of_neighbors(NeighborhoodSearch const& nsearch)
{
	unsigned long res = 0;
	auto const& d = nsearch.point_set(0);
	for (int i = 0; i < d.n_points(); ++i)
	{
		res += static_cast<unsigned long>(d.n_neighbors(i));
	}
	return static_cast<Real>(res) / static_cast<Real>(d.n_points());
}

Real
compute_average_distance(NeighborhoodSearch const& nsearch)
{
	unsigned long long res = 0;
	auto const& d = nsearch.point_set(0);
	unsigned long long count = 0;
	for (int i = 0; i < d.n_points(); ++i)
	{
		std::size_t nn = d.n_neighbors(i);
		for (int j = 0; j < nn; ++j)
		{
			CompactNSearch::PointID const& k = d.neighbor(i, j);
			res += std::abs(i - static_cast<int>(k.point_id));
			count++;
		}
	}
	return static_cast<Real>(res) / static_cast<Real>(count);
}

std::vector<std::vector<unsigned int>>
brute_force_search()
{
	std::vector<std::vector<unsigned int>> brute_force_neighbors(positions.size());
	for (int i = 0; i < positions.size(); ++i)
	{
		std::vector<unsigned int>& neighbors = brute_force_neighbors[i];
		for (int j = 0; j < positions.size(); ++j)
		{
			if (i == j)
				continue;
			std::array<Real, 3> const& xa = positions[i];
			std::array<Real, 3> const& xb = positions[j];
			Real l2 =
					(xa[0] - xb[0])*(xa[0] - xb[0]) +
					(xa[1] - xb[1])*(xa[1] - xb[1]) +
					(xa[2] - xb[2])*(xa[2] - xb[2]);
			if (l2 <= radius * radius)
			{
				neighbors.push_back(j);
			}
		}
	}
	return std::move(brute_force_neighbors);
}

void
compare_with_bruteforce_search(NeighborhoodSearch const& nsearch)
{
	auto brute_force_neighbors = brute_force_search();
	PointSet const& d0 = nsearch.point_set(0);
	for (int i = 0; i < N; ++i)
	{
		auto const& bfn = brute_force_neighbors[i];
		if (bfn.size() != d0.n_neighbors(i))
		{
			std::cerr << "ERROR: Not the same number of neighbors." << std::endl;
		}
		for (int j = 0; j < d0.n_neighbors(i); ++j)
		{
			if (std::find(bfn.begin(), bfn.end(), d0.neighbor(i, j).point_id) == bfn.end())
			{
				std::cerr << "ERROR: Neighbor not found in brute force list." << std::endl;
			}
		}
	}
}

std::array<Real, 3>
enright_velocity_field(std::array<Real, 3> const& x)
{
	Real sin_pi_x_2 = std::sin(M_PI * x[0]);
	Real sin_pi_y_2 = std::sin(M_PI * x[1]);
	Real sin_pi_z_2 = std::sin(M_PI * x[2]);
	sin_pi_x_2 *= sin_pi_x_2;
	sin_pi_y_2 *= sin_pi_y_2;
	sin_pi_z_2 *= sin_pi_z_2;

	Real sin_2_pi_x = std::sin(2.0 * M_PI * x[0]);
	Real sin_2_pi_y = std::sin(2.0 * M_PI * x[1]);
	Real sin_2_pi_z = std::sin(2.0 * M_PI * x[2]);
	return {{
			2.0 * sin_pi_x_2 * sin_2_pi_y * sin_2_pi_z,
			-sin_2_pi_x * sin_pi_y_2 * sin_2_pi_z,
			-sin_2_pi_x * sin_2_pi_y * sin_pi_z_2}};

}

void
advect()
{
#ifdef _MSC_VER
	concurrency::parallel_for_each
#else
	__gnu_parallel::for_each
#endif
	(positions.begin(), positions.end(), [&](std::array<Real, 3>& x)
	{
		std::array<Real, 3> v = enright_velocity_field(x);
		x[0] += 0.005 * v[0];
		x[1] += 0.005 * v[1];
		x[2] += 0.005 * v[1];
	}
	);
}

int main(int argc, char* argv[])
{
	Real min_x = std::numeric_limits<Real>::max();
	Real max_x = std::numeric_limits<Real>::min();
	positions.reserve(N * N * N);
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < N; ++j)
		{
			for (unsigned int k = 0; k < N; ++k)
			{
				std::array<Real, 3> x = {{
						r_omega * (2.0 * static_cast<Real>(i) / static_cast<Real>(N-1)-1.0),
						r_omega * (2.0 * static_cast<Real>(j) / static_cast<Real>(N-1)-1.0),
						r_omega * (2.0 * static_cast<Real>(k) / static_cast<Real>(N-1)-1.0)}};

				Real l2  = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
				if (l2 < r_omega2)
				{
					x[0] += 0.35;
					x[1] += 0.35;
					x[2] += 0.35;
					positions.push_back(x);
					if (min_x > x[0])
					{
						min_x = x[0];
					}
					if (max_x < x[0])
					{
						max_x = x[0];
					}
				}
			}
		}
	}
	std::random_shuffle(positions.begin(), positions.end());

	NeighborhoodSearch nsearch(radius, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.find_neighbors();

	std::cout << "#Points                                = " << positions.size() << std::endl;
	std::cout << "Search radius                          = " << radius << std::endl;
	std::cout << "Min x                                  = " << min_x << std::endl;
	std::cout << "Max x                                  = " << max_x << std::endl;
	std::cout << "Average number of neighbors            = " << compute_average_number_of_neighbors(nsearch) << std::endl;
	std::cout << "Average index distance prior to z-sort = " << compute_average_distance(nsearch) << std::endl;

	nsearch.z_sort();
	for (auto const& d : nsearch.point_sets())
	{
		d.sort_field(positions.data());
	}
	nsearch.find_neighbors();
	//compare_with_bruteforce_search(nsearch);

	std::cout << "Average index distance after z-sort    = " << compute_average_distance(nsearch) << std::endl;

	std::cout << "Moving points:" << std::endl;
	for (int i = 0; i < N_enright_steps; ++i)
	{
		std::cout << "Enright step " << i << ". ";
		advect();
		auto t0 = std::chrono::high_resolution_clock::now();
		nsearch.find_neighbors();
		std::cout << "Neighborhood search took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << "ms" << std::endl;
		//compare_with_bruteforce_search(nsearch);
	}

	return 0;
}
