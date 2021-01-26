#include "stdafx.h"
#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include "iostream"
#include "iomanip"
#include <math.h>
using namespace std;
const double h1 = 0.1, h2 = 0.01;

double funcz(double x, double u, double z)
{
	return (2 * exp(-x)*cos(x)-2*z-2*u)
}
double func(double x)
{
	return exp(-x)*(cos(x) + x*sin(x) + sin(x));
}
int _tmain(int argc, _TCHAR* argv[])
{
	FILE *fout;
	errno_t err;
	err = fopen_s(&fout, "Cm_4.txt", "wt");
	int i;
	double x, u = 1, z = 0, u1, z1, uu = 1, zz = 0;
	fprintf(fout, "%0.5f;%0.5f; %0.5f;\n", u, uu, func(0));
	for (i = 1; i <= 10; i++)
	{
		u = 1; z = 0;
		for (x = 0; x < i*0.1 - 0.001; x = x + h1)
		{
			u1 = u + h1*z + h1*h1*funcz(x, u, z) / 2;
			z1 = z + h1*(funcz(x, u, z) + funcz(x + h1, u1, z + h1*funcz(x, u, z))) / 2;
			u = u1; z = z1;
		}
		uu = 1; zz = 0;
		for (x = 0; x < i*0.1 - 0.0001; x = x + h2)
		{
			u1 = uu + h2*zz + h2*h2*funcz(x, uu, zz) / 2;
			z1 = zz + h2*(funcz(x, uu, zz) + funcz(x + h2, u1, zz + h2*funcz(x, uu, zz))) / 2;
			uu = u1; zz = z1;
		}

		fprintf(fout, "%0.5f;%0.5f; %0.5f;\n", u, uu, func(i*0.1));
	}
	system("pause");
	return 0;
}
