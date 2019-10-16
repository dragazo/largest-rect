#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <utility>
#include <array>
#include <sstream>
#include <chrono>

#include "rang/include/rang.hpp"

std::mt19937 rng{ std::random_device{}() }; // create a high-quality pseudo random number generator

unsigned int rand_bool_chance = 2;
bool rand_bool() { return rng() % rand_bool_chance; }

typedef std::vector<std::vector<bool>>      matrix_t;
typedef std::pair<std::size_t, std::size_t> vec2_t;
typedef std::pair<vec2_t, vec2_t>           rectangle_t;

// returns a rows x cols matrix of bools
matrix_t rand_matrix(std::size_t rows, std::size_t cols)
{
	matrix_t res(rows);
	for (auto &row : res) std::generate_n(std::back_inserter(row), cols, rand_bool);
	return res;
}
// prints the matrix to stdout
void print_solution(const matrix_t &matrix, const rectangle_t &solution, bool verbose)
{
	std::cout << "solution: pos=(" << solution.first.first << ", " << solution.first.second
		<< ") size=(" << solution.second.first << ", " << solution.second.second << ")=" << (solution.second.first * solution.second.second) << '\n';
	
	if (verbose)
	{
		std::cout << rang::style::bold << rang::style::reversed;

		for (std::size_t i = 0; i < matrix.size(); ++i)
		{
			for (std::size_t j = 0; j < matrix[i].size(); ++j)
			{
				if (i >= solution.first.first && i < solution.first.first + solution.second.first && j >= solution.first.second && j < solution.first.second + solution.second.second) std::cout << rang::fg::blue;
				else if (matrix[i][j]) std::cout << rang::fg::green;
				else std::cout << rang::fg::red;

				std::cout << matrix[i][j] << ' ';
			}
			std::cout << '\n';
		}

		std::cout << rang::style::reset;
	}
}

// finds the largest rectangle in the given histogram (defined as an array of heights).
// returns an array of values: [0]=left index of rectangle, [1]=width of rectangle, [2]=height of rectangle
std::array<std::size_t, 3> find_max_histogram_area(const std::vector<std::size_t> &histogram)
{
	std::vector<std::size_t> stack; // we preform our logic from left to right by storing indexes from histogram array (initially empty)
	std::size_t i = 0;              // working index
	std::size_t best_area = 0;      // keep track of the best area we've managed to find
	std::size_t best_left = 0;      // keep track of the left starting position of the largest rectangle
	std::size_t best_width = 0;     // keep track of the width of the biggest rectangle
	std::size_t best_height = 0;    // keep track of the height of the biggest rectangle

	auto pop_one = [&] // helper function - pops one item from the stack and performs max logic
	{
		std::size_t top = stack.back(); // get the top item on the stack and pop it off
		stack.pop_back();

		// compute area using top bar we just popped off as the smallest
		// if stack is not empty, we go up to and including the new top of stack index, otherwise we use everything to the left (popped current minimum height)
		std::size_t area = histogram[top] * (stack.empty() ? i : i - 1 - stack.back());

		// update best area info
		if (area > best_area)
		{
			best_area = area;
			best_left = stack.empty() ? 0 : 1 + stack.back();
			best_width = stack.empty() ? i : i - 1 - stack.back();
			best_height = histogram[top];
		}
	};

	// for each column of the histogram
	while (i < histogram.size())
	{
		// if the stack is currently empty or if this value is at least as large as the previous, add it to the stack and move to the next column
		if (stack.empty() || histogram[stack.back()] <= histogram[i]) stack.push_back(i++);
		// otherwise stack is not empty and this value is smaller than stack top - perform pop one logic
		else pop_one();
	}
	// then run pop one logic until the stack is empty
	while (!stack.empty()) pop_one();

	return { best_left, best_width, best_height };
}

rectangle_t find_largest_rect_dynamic(const matrix_t &matrix)
{
	std::vector<std::size_t> histogram(matrix[0].size(), 0); // create a 1D array to represent histogram info (fill with zeros initially)
	rectangle_t best = { { 0, 0 }, { 0, 0 } }; // keep track of the best solution
	std::size_t best_area = 0;                 // as well as the actual area to avoid a bunch of multiplications

	// for each row in the matrix
	for (std::size_t i = 0; i < matrix.size(); ++i)
	{
		// update the histogram info (1 will add 1 to height, 0 will reset height to 0)
		for (std::size_t j = 0; j < matrix[i].size(); ++j)
			if (matrix[i][j]) histogram[j] += 1; else histogram[j] = 0;

		auto res = find_max_histogram_area(histogram); // get the largest area from the histogram
		std::size_t area = res[1] * res[2]; // compute the area of the result from the histogram

		// update best rectangle info
		if (area > best_area)
		{
			best_area = area;
			best = { { i - res[2] + 1, res[0] }, { res[2], res[1] } };
		}
	}

	return best; // return the best solution we managed to find
}

// returns true if the given solution (rectangle) is a valid solution for the given matrix.
// returning true means the rectangle is within bounds and contains only 1's.
bool is_valid(const matrix_t &matrix, const rectangle_t &solution)
{
	// if the solution is out of bounds of the matrix, it's invalid
	if (solution.first.first + solution.second.first > matrix.size()) return false;
	if (solution.first.second + solution.second.second > matrix[0].size()) return false;

	// if the solution contains a 0 it's invalid
	for (std::size_t i = 0; i < solution.second.first; ++i)
		for (std::size_t j = 0; j < solution.second.second; ++j)
			if (!matrix[solution.first.first + i][solution.first.second + j]) return false;

	// otherwise it's valid
	return true;
}

// finds the largest rectangle of 1's in the matrix using a branch and bound algorithm.
// due to its simplicity, the recursion is unfolded to be a simple loop.
// because if this, it might actually just be a brute force algorithm with pruning - depends on technical definitions of both.
rectangle_t find_largest_rect_branch_bound(const matrix_t &matrix)
{
	rectangle_t solution;                           // the working solution
	rectangle_t best_solution = { {0, 0}, {0, 0} }; // the best solution we've come up with so far
	std::size_t best_size = 0;                      // keep track of the size of the best solution separately to avoid a bunch of multiplies

	// for each starting row position
	for (solution.first.first = 0; solution.first.first < matrix.size(); ++solution.first.first)
	{
		// and for each starting col position
		for (solution.first.second = 0; solution.first.second < matrix[0].size(); ++solution.first.second)
		{
			// and over all vertical sizes
			for (solution.second.first = 1; solution.first.first + solution.second.first-1 < matrix.size() && matrix[solution.first.first + solution.second.first-1][solution.first.second]; ++solution.second.first)
			{
				// and over all horizontal sizes
				for (solution.second.second = 1; is_valid(matrix, solution); ++solution.second.second)
				{
					std::size_t size = solution.second.first * solution.second.second;
					if (size > best_size)
					{
						best_solution = solution;
						best_size = size;
					}
				}
			}
		}
	}

	return best_solution;
}

// parses a value of type T from the string and returns true iff it was successful
template<typename T>
bool parse(T &dest, const char *str)
{
	std::istringstream istr{ str };
	istr >> dest;
	return istr && istr.get() == EOF;
}
// as parse() except that failure triggers an error message and immediate program termination
template<typename T, int errorcode = 200>
T volatile_parse(const char *str)
{
	T temp;
	if (!parse(temp, str)) { std::cerr << "failed to parse " << str << " as " << typeid(T).name() << '\n'; std::exit(errorcode); }
	return temp;
}

int main(int argc, const char *const argv[])
{
	using namespace std::chrono;

	if (argc < 4)
	{
		std::cerr << "usage: " << argv[0] << " [zero chance] [rows] [cols]\n";
		return 1;
	}

	int zero_chance = volatile_parse<int>(argv[1]);
	if (zero_chance <= 0) { std::cerr << "attempt to set zero chance to " << zero_chance << " (must be positive)\n"; return 2; }
	rand_bool_chance = zero_chance;

	int rows = volatile_parse<int>(argv[2]);
	int cols = volatile_parse<int>(argv[3]);
	if (rows <= 0 || cols <= 0) { std::cerr << "attempt to use empty matrix\n"; return 2; }

	bool verbose = false; // default verbose to false
	
	// process extra options
	for (int i = 4; i < argc; ++i)
	{
		if (strcmp(argv[i], "-v") == 0) verbose = true;
		else { std::cerr << "unknown command line option: " << argv[i] << '\n'; return 3; }
	}

	auto matrix = rand_matrix(rows, cols);

	struct
	{
		const char *name;                       // name of the algorithm
		rectangle_t(*perform)(const matrix_t&); // function to perform the algorithm and get the result
	} algorithms[] = 
	{
		{ "dynamic", find_largest_rect_dynamic },
		{ "branch and bound", find_largest_rect_branch_bound },
	};

	for (const auto &algorithm : algorithms)
	{
		// perform the search using this algorithm and time how long it took to complete (in nanoseconds)
		auto start = high_resolution_clock::now();
		auto solution = algorithm.perform(matrix);
		auto stop = high_resolution_clock::now();
		auto timeus = duration_cast<microseconds>(stop - start).count();

		// print the results for this algorithm's solution
		std::cout << algorithm.name << ": " << ((double)timeus / 1000000) << " sec (" << timeus << " us)\n";
		print_solution(matrix, solution, verbose);
		std::cout << std::endl;
	}

	return 0;
}