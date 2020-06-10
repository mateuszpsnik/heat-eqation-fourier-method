/*
	Mateusz Psnik
	170636
*/

#include <iostream>
#include <cmath>
#include <fstream>

const double pi = 3.14159;

double integral(double lambda, double lower_bound = 0, double upper_bound = 2 * pi)
{
	double sum = 0.0;
	const double step = 0.2;

	for (double x = step / 2; x < upper_bound; x += step)
	{
		sum += step * std::sin(x) * std::sin(lambda * x);
	}

	return sum;
}

double g(double x, double t)
{
	return 2 * std::exp(-t) * std::sin(x);
}

double integral_g_n(double lambda, double t_i, double L = 2 * pi)
{
	double sum = 0.0;
	const double step = 0.05;
	const int NX = (int)(L / step);

	for (size_t i = 0; i < NX; i++)
	{
		sum += g(i * step, t_i) * std::sin(lambda * i * step);
	}

	return (double)(2.0 / NX) * sum;
}

double integral_t_n(double lambda, double t, double L = 2 * pi, double kappa = 1)
{
	double sum = 0.0;
	const double tau = 0.05;
	const int M = (int)(t / tau);

	for (size_t i = 0; i < M; i++)
	{
		double a = std::exp(kappa * lambda * lambda * (tau * i - t));
		double b = integral_g_n(lambda, tau * i);

		sum += std::exp(kappa * lambda * lambda * (tau * i - t))
			* integral_g_n(lambda, tau * i) * tau;
	}

	return sum;
}

double solve_nonhom_eq(double x, double t, int n, double kappa = 1.0, double L = 2 * pi)
{
	double sum = 0.0;
	double lambda = 0.0;
	
	for (size_t i = 0; i <= n; i++) 
	{
		lambda = i * pi / L;

		sum += integral_t_n(lambda, t) * std::sin(lambda * x);
	}

	return sum;
}

double solve_hom_eq(double x, double t, int n, double kappa = 1.0, double L = 10.0)
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

void solve_write(std::string filename, double t, int n = 100, double kappa = 1.0, double L = 2 * pi)
{
	using std::endl;
	std::ofstream file{ filename };

	file << "t = " << t << endl;
	file << "n = " << n << endl;

	double step = 0.05;

	for (double x = 0.0; x <= L; x += step)
	{
		double solution = solve_hom_eq(x, t, n) + solve_nonhom_eq(x, t, n);
		file << x << "," << solution << endl;
	}
}

int main()
{
	std::cout << "Writing the result into text files might take some time (about 1-2 minutes)\n";

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
}