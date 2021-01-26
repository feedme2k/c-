// CM_AB_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <tchar.h>
#include "iostream"
#include "iomanip"
#include <math.h>
using namespace std;
double eps = 1.0e-6;

double Fx1(double x, double y)
{
	return sin(x + y) - (1.1 * x) + 0.2;
}
double Fx2(double x, double y)
{
	return pow(x, 2) + pow(y, 2) - 1;
}
double dFx1dx(double x, double y)
{
	return cos(x + y) - 1.1;
}
double dFx1dy(double x, double y)
{
	return cos(x + y);
}
double dFx2dx(double x, double y)
{
	return 2 * x;
}
double dFx2dy(double x, double y)
{
	return 2 * y;
}
double det(double a[4])
{
	return a[0] * a[3] - a[1] * a[2];
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << setprecision(7);
	int i = 0;
	double x[2] = { 0.5, 0.5 }, x1[2], W[4], W1[4], K;
	cout << "x[0] = " << x[0] << " " << "x[1] = " << x[1] << " " << Fx1(x[0], x[1]) << " " << Fx2(x[0], x[1]);
	system("pause");
	W[0] = dFx1dx(x[0], x[1]);
	W[1] = dFx1dy(x[0], x[1]);
	W[2] = dFx2dx(x[0], x[1]);
	W[3] = dFx2dy(x[0], x[1]);
	K = det(W);
	W1[0] = W[3] / K;
	W1[1] = -W[1] / K;
	W1[2] = -W[2] / K;
	W1[3] = W[0] / K;

	x1[0] = x[0] - W1[0] * Fx1(x[0], x[1]) - W1[1] * Fx2(x[0], x[1]);
	x1[1] = x[1] - W1[2] * Fx1(x[0], x[1]) - W1[3] * Fx2(x[0], x[1]);

	while (abs(x1[0] - x[0]) > eps && abs(x1[1] - x[1]) > eps && abs(Fx1(x1[0], x1[1])) > eps && abs(Fx2(x1[0], x1[1])) > eps)
	{
		x[0] = x1[0];
		x[1] = x1[1];
		W[0] = dFx1dx(x[0], x[1]);
		W[1] = dFx1dy(x[0], x[1]);
		W[2] = dFx2dx(x[0], x[1]);
		W[3] = dFx2dy(x[0], x[1]);
		K = det(W);
		W1[0] = W[3] / K;
		W1[1] = -W[1] / K;
		W1[2] = -W[2] / K;
		W1[3] = W[0] / K;
		x1[0] = x[0] - W1[0] * Fx1(x[0], x[1]) - W1[1] * Fx2(x[0], x[1]);
		x1[1] = x[1] - W1[2] * Fx1(x[0], x[1]) - W1[3] * Fx2(x[0], x[1]);
		i++;
	}
	cout << "x=" << x1[0] << " y=" << x1[1] << endl << " f1(x,y)=" << Fx1(x1[0], x1[1]) << " f2(x,y)=" << Fx2(x1[0], x1[1]) << endl;
	cout << "Количество итераций i=" << i << endl;
	system("pause");
	return 0;
}
