#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>

std::mt19937 rng{ std::random_device{}() }; // create a high-quality pseudo random number generator

// returns a random bool value (50/50)
bool rand_bool() { return rng() & 1; }

typedef std::vector<std::vector<bool>> matrix_t;

// returns a rows x cols matrix of bools
auto rand_matrix(std::size_t rows, std::size_t cols)
{
	matrix_t res(rows);
	for (auto &row : res) std::generate_n(std::back_inserter(row), cols, rand_bool);
	return res;
}
// prints the matrix to stdout
void print(const matrix_t &matrix)
{
	for (const auto &row : matrix)
	{
		for (auto i : row) std::cout << i << ' ';
		std::cout << '\n';
	}
}

int main()
{
	auto matrix = rand_matrix(25, 50);

	print(matrix);

	std::cin.get();
	return 0;
}