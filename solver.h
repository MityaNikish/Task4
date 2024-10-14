#pragma once
#include <vector>

struct Values
{
	//�������� � �����e
	double p;
	//��������� � �����e
	double ro;
	//�������� � �����e
	double u;
};

class Solver
{
	//�������� �� ����� � �����
	double p0_;
	//��������� �� ����� � �����
	double ro0_;
	//������� ������� ������
	std::vector<double> A_;

	//����������� �������
	double x_star_;
	//������� ������������ �������
	double A_star_;

	//���������� ������ ����������
	double x_s_;

	double gamma_ = 2;

public:
	Solver(const double p0, const double ro0, const std::vector<double>& A, const double x_s = 1);
	std::vector<Values> solving();

private:
	std::pair<double, double> searchMin();
};