#include <iostream>
#include <cmath>

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

double solve_eq(double x, double t, int n, double kappa = 1, double L = 10)
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

int main()
{
	


}