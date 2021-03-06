/*
	Mateusz Psnik
	170636
*/

#include <iostream>
#include <cmath>
#include <fstream>

double integral(double lambda, double lower_bound = 0, double upper_bound = 10)
{
	double sum = 0.0;
	const double step = 0.2;

	for (double x = step / 2; x < upper_bound; x += step)
	{
		sum += step * std::exp(-(x - 5) * (x - 5)) 
			* std::sin(lambda * x);
	}

	return sum;
}

double solve_eq(double x, double t, int n, double kappa = 1.0, double L = 10.0)
{
	double sum = 0.0;
	double lambda = 0.0;
	const double pi = 3.14159;

	for (size_t i = 0; i <= n; i++)
	{
		lambda = i * pi / L;

		sum += 2 / L * integral(lambda) * std::exp(-kappa * lambda * t)
			* std::sin(lambda * x);
	}

	return sum;
}

void solve_write(std::string filename, double t, int n = 100, double kappa = 1.0, double L = 10.0)
{
	using std::endl;
	std::ofstream file{ filename };

	file << "t = " << t << endl;
	file << "n = " << n << endl;

	double step = 0.05;

	for (double x = 0.0; x <= L; x += step)
	{
		file << x << "," << solve_eq(x, t, n) << endl;
	}
}

int main()
{
	solve_write("05_5.txt", .5, 5);
	solve_write("05_20.txt", .5, 20);
	solve_write("05_80.txt", .5, 80);

	solve_write("01.txt", 0.1);
	solve_write("05.txt", 0.5);
	solve_write("1.txt", 1);
	solve_write("2.txt", 2);
	solve_write("4.txt", 4);
	solve_write("6.txt", 6);
	solve_write("8.txt", 8);
	solve_write("10.txt", 10);
}