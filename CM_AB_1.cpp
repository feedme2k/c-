// CM_AB_1.cpp : Defines the entry point for the console application.
//

#include <stdafx.h>
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <math.h>
#include <conio.h>
#include <limits.h>
#include <iomanip>
using namespace std;

double f(double x)
{
	return 1 - x + sin(x) - log(1 + x);
}

void main(){
	setlocale(LC_CTYPE, "Russian");
	double a, b, eps, L, R, fL, fR, m, fm;
	int n;
	cout << "\n ������� �������� eps \n" << endl;
	cin >> eps;
	cout << "\n ����� ���� \n" << endl;
	cin >> a;
	cout << "\n ������ ���� \n" << endl;
	cin >> b;

	cout << setprecision(8);
	//����� ������� ������� ������� 

	L = a; R = b;
	fL = f(L);
	fR = f(R);
	n = 0; // ���������� �������� 
	while (fabs(R - L)> eps) {
		m = (R + L) / 2;
		n++;
		fm = f(m);
		if (fm == 0)
			break;
		if (fL*fm <= 0){
			R = m; fR = fm;
		}
		else {
			L = m; fL = fm;
		}
	}
	m = (L + R) / 2;
	cout << "\n ��������� �� ������ ������� ������� \n";
	cout << " ����� �������� \n" << " " << n << " " << "���������" << " " << m << "\n";

	// ����� ���� 
	L = a; R = b; n = 0;
	fL = f(L); fR = f(R);
	while (fabs(R - L)> eps) {
		m = L - fL*(R - L) / (fR - fL);
		n++;
		fm = f(m);
		if (fm == 0)
			break;
		if (fL*fm <= 0) {
			R = m; fR = fm;
		}
		else {
			L = m; fL = fm;
		}
	}
	m = (L + R) / 2.0;
	cout << "\n ����� ��������=" << " " << n << " " << "\n ���������=" << " " << m << "\n";
	getch;
	system("pause");
}

