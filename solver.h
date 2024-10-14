#pragma once
#include <vector>

struct Values
{
	//Давление в каналe
	double p;
	//Плотность в каналe
	double ro;
	//Скорость в каналe
	double u;
};

class Solver
{
	//Давление на входе в канал
	double p0_;
	//Плотность на входе в канал
	double ro0_;
	//Площадь сечения канала
	std::vector<double> A_;

	//Критическое сечение
	double x_star_;
	//Площадь критического сечения
	double A_star_;

	//Координаты скачка уплотнения
	double x_s_;

	double gamma_ = 2;

public:
	Solver(const double p0, const double ro0, const std::vector<double>& A, const double x_s = 1);
	std::vector<Values> solving();

private:
	std::pair<double, double> searchMin();
};