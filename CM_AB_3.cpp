// CM_AB_3.cpp : Defines the entry point for the console application.
//

#include <stdafx.h>
#include <iostream>
#include <cmath>
#include <conio.h>
#include <iomanip>
#include <math.h>
#define EPS 1e-7
#define M_PI 3.1415926535897932384626433832795
using namespace std;

double f(double x)
{
	return pow(tan(x), 2) + pow(1 / tan(x), 2) - (2 * x);
}

double fth(double x)
{
	return tan(x)-(1/tan(x))-(2*x);
}

double simpson(double a, double b, long N)
{
	double I, I2 = 0, I4 = 0;
	double step = fabs(b - a) / N;
	I = f(a) + f(b);
	for (double k = a + step; k < b; k += 2 * step)
		I4 += 4 * f(k);
	for (double k = a + 2 * step; k < b; k += 2 * step)
		I2 += 2 * f(k);
	I += I4 + I2;
	return (I*(step / 3));

}
double Integral(double a, double b){

	int N = 16;
	double In, I2n, D;
	In = simpson(a, b, N);
	do{
		N *= 2;
		I2n = simpson(a, b, N);
		D = (I2n - In);
		In = I2n;
	} while (fabs(D) >= EPS);
	return I2n + D;
}


int main()
{
	double x1, x2, Fx = 0, FxTh, A = M_PI/36, B = (17*M_PI)/36, h =M_PI/36;
	cout << setprecision(7);
	cout << setiosflags(ios_base::fixed);
	Fx = Integral(A, A);
	FxTh = fth(A);
	cout << "[X]:" << setw(3) << A << " [F(X)]:" << setw(7) << Fx << " [Ft(x)]:" << setw(7) << FxTh << endl;
	x2 = A;
	for (int i = 0; i*h <= B - A; i++)
	{
		x1 = x2;
		x2 = x1 + h;
		Fx += Integral(x1, x2);
		FxTh = fth(x2);
		cout << "[X]:" << setw(3) << x2 << " [F(X)]:" << setw(7) << Fx << " [Ft(x)]:" << setw(7) << FxTh << endl;
	}
	_getch();

}
