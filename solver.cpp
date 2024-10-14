#include "solver.h"

#include <utility>
#include <functional>
#include <cfloat>

namespace 
{
	double halfCutMethod(const std::function<double(double)>& func, const double left_boundary, const double right_boundary, const double epsilon, const size_t quantity_iteration)
	{
		double a = left_boundary;
		double b = right_boundary;
		double c = (a + b) / 2;

		for (size_t i = 0; i < quantity_iteration; ++i)
		{
			double f_a = func(a);
			double f_b = func(b);

			c = (a + b) / 2;
			double f_c = func(c);

			a = (f_a * f_c) > 0 ? c : a;
			b = (f_b * f_c) > 0 ? c : b;

			if (abs(f_c) < epsilon)
			{
				break;
			}
		}
		return c;
	}
}

Solver::Solver(const double p0, const double ro0, const std::vector<double>& A, const double x_s) : p0_(p0), ro0_(ro0), A_(A), x_s_(x_s)
{
	const std::pair<double, double> min_values = searchMin();
	x_star_ = min_values.first;
	A_star_ = min_values.second;
}

std::vector<Values> Solver::solving()
{
	const double gamma = gamma_;
	const double A_star = A_star_;
	const double A_0 = A_[0];

	std::function<double(double)> func = [gamma, A_0, A_star, M_star = 1] (double M)
	{
		const double Q = (1 + (gamma - 1) / 2 * M_star) / (1 + (gamma - 1) / 2 * M);
		return M_star / M * pow(Q, (gamma + 1) / (1 - gamma) / 2) * A_star - A_0;
	};

	const double M_0 = halfCutMethod(func, 0.0001, 1, 0.0001, 100000);
	const double u0 = pow(gamma * p0_ / ro0_, 0.5);



	const double h = 1 / A_.size();
	double M;
	double Q;

	std::vector<Values> data(A_.size());

	for (size_t i = 0; i < A_.size(); ++i)
	{
		const double A_i = A_[i];
		func = [gamma, A_i, A_0, M_0](double M)
			{
				double Q = (1 + (gamma - 1) / 2 * M_0) / (1 + (gamma - 1) / 2 * M);
				return M_0 / M * pow(Q, (gamma + 1) / (1 - gamma) / 2) * A_0 - A_i;
			};

		if (i * h <= x_star_)
		{
			M = halfCutMethod(func, 0.0001, 1, 0.0001, 100000);
		}
		else
		{
			M = halfCutMethod(func, 1, 10, 0.0001, 100000);
		}

		Q = (1 + (gamma - 1) / 2 * M_0) / (1 + (gamma - 1) / 2 * M);

		data[i].ro = ro0_ * pow(Q, 1 / (gamma - 1));
		data[i].p = p0_ * pow(Q, gamma / (gamma - 1));
		data[i].u = u0 * pow(Q, 0.5) * M / M_0;
	}

	return data;

}

std::pair<double, double> Solver::searchMin()
{
	double A_min = DBL_MAX;
	size_t index_min = 0;

	for (size_t i = 0; i < A_.size(); ++i)
	{
		if (A_[i] < A_min)
		{
			index_min = i;
			A_min = A_[i];
		}
	}
	double x_min = index_min / A_.size();

	return { x_min, A_min };
}