// Vetvei i granic.cpp : Defines the entry point for the console application.
//

#include"stdafx.h"
#include<iostream>
#include<conio.h>
#include<math.h>
#include<iomanip>
#include<limits.h>
using namespace std;
const int def = 5;

const bool debug = true;

const double M = 2000;
//Branch_and_Bound
struct BranchAndBound{
	double**Mx;
	int n, m;
	int f;
	BranchAndBound *l, *r;
};

bool isMaximumOptimal(double ** matr, int n, int m){
	for (int i = 0; i<m - 1; i++)
	if (matr[n - 1][i] < 0)
		return 0;
	return 1;
}
bool isMinimumOptimal(double ** matr, int n, int m){
	for (int i = 0; i<m - 1; i++)
	if (matr[n - 1][i] > 0)
		return 0;
	return 1;
}
bool isSolutionsExistence(double ** matr, int n, int m, int numb){
	for (int i = 0; i<n - 1; i++)
	if (matr[i][numb] > 0)
		return 0;
	return 1;
}
void divisionLine(double**matr, int n, int m, int column, int line){
	double a = matr[line][column];
	for (int i = 0; i<m; i++)
		matr[line][i] /= a;
}
void subtractionLineWithCoefficient(double**matr, int n, int m, int b, int row, int c){
	double a = matr[row][b];
	for (int i = 0; i<m; i++)
		matr[row][i] -= a*matr[c][i];
}
void add(double**matr, int n, int m, int numb, double w, int numb1){
	for (int i = 0; i<m; i++)
		matr[numb][i] += w*matr[numb1][i];
}
bool isArtificialVariableInBasis(double*res, int n, int m){
	int i;
	for (i = m - n; i < m - 1 && res[i] == 0; i++);
	if (i == m - 1) return 0;
	else return 1;
}

double fractionalPart(double x){
	if (x>0 && x - floor(x) <= 0.99999)
		return x - floor(x);
	else if (x<0)
		return fabs(ceil(x)) + 1 - fabs(x);
	else
		return 0;
}
double** inputMatrix(int n, int m)
{
	double**matr = new double*[n];
	for (int i = 0; i<n; i++)
	{
		matr[i] = new double[m];
		for (int j = 0; j<m; j++)
			matr[i][j] = 0;
	}
	return matr;
}

double* inputArray(int m)
{
	double*res = new double[m - 1];
	for (int i = 0; i<m - 1; i++)
		res[i] = 0;
	return res;
}

void equateMatrix(double**m1, double**m2, int n, int m)
{
	for (int i = 0; i<n; i++)
	for (int j = 0; j<m; j++)
		m2[i][j] = m1[i][j];
}

double** deployNewMatrix(double**matr, int n, int m, double numb, int val, bool flag, bool var)
{
	double**mx = inputMatrix(n + 1, m + n + 1);
	for (int i = 0; i<n - 1; i++)
	for (int j = 0; j<m - 1; j++)
		mx[i][j] = matr[i][j];
	for (int i = 0; i<n - 1; i++)
		mx[i][m + n] = matr[i][m - 1];
	for (int j = 0; j<m - 1; j++)
	{
		mx[n][j] = matr[n - 1][j];
	}
	mx[n - 1][val] = 1;
	if (flag) mx[n - 1][m - 1] = 1;
	else mx[n - 1][m - 1] = -1;
	if (numb == 0) mx[n - 1][m - 1] = 0;
	mx[n][m + n] = matr[n - 1][m - 1];
	mx[n - 1][m + n] = numb;
	for (int i = 0; i<n; i++)
		mx[i][m + i] = 1;
	if (var)
	for (int j = m; j<m + n; j++)
		mx[n][j] = M;
	else
	for (int j = m; j<m + n; j++)
		mx[n][j] = -M;
	return mx;
}

bool isResulfInteger(double*res, int m, int*y, int k)
{
	bool f = 1;
	for (int j = 0; j<k; j++)
	if (fractionalPart(res[y[j] - 1])<pow((double)10, -5))
		f = 1;
	else {
		f = 0;
		break;
	}
	if (f)
		return 1;
	else
		return 0;
}


double * prepareResult(double**Matrix, int n, int m)
{
	double * res = inputArray(m - 1);
	int i, j, l, k, ind;
	for (j = 0; j<m - 1; j++)
	{
		l = 0; k = 0;
		for (i = 0; i<n; i++)
		if ((int)fabs(Matrix[i][j]) == 0)
			l++;
		else if ((int)Matrix[i][j] == 1) {
			k++;
			ind = i;
		}
		if (l == n - 1 && k == 1)
			res[j] = Matrix[ind][m - 1];
		else
			res[j] = 0;
	}
	return res;
}
bool Simplex_Method(double**Matrix, int n, int m, bool var)
{
	double max = INT_MIN, min = INT_MAX, min_otn = INT_MAX, p;
	int index_v_basis, index_iz_basis, i;

	if (var){
		for (i = 0; i<m - 1; i++)
		if (Matrix[n - 1][i]<min) { min = Matrix[n - 1][i]; index_v_basis = i; }
	}
	else {
		for (i = 0; i<m - 1; i++)
		if (Matrix[n - 1][i]>max) { max = Matrix[n - 1][i]; index_v_basis = i; }
	}
	if (!isSolutionsExistence(Matrix, n, m, index_v_basis)){
		for (i = 0; i<n - 1; i++)
		{
			if (Matrix[i][index_v_basis]>0)
			{
				p = Matrix[i][m - 1] / Matrix[i][index_v_basis];
				if (p<min_otn) { min_otn = p; index_iz_basis = i; }
			}
		}
		divisionLine(Matrix, n, m, index_v_basis, index_iz_basis);
		for (i = 0; i<n; i++)
		if (i != index_iz_basis) subtractionLineWithCoefficient(Matrix, n, m, index_v_basis, i, index_iz_basis);
		return 1;
	}
	else return 0;
}

void printSimplexTable(double**Matrix, int n, int m){
	for (int k = 0; k<n; k++){
		for (int j = 0; j<m; j++)
			cout << setw(3) << setprecision(3) << " " << Matrix[k][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void printResulfLine(double**Matrix, double*res, int n, int m, bool f){
	cout << "Decision: ";
	for (int i = 0; i<m - 1; i++)
		cout << setprecision(4) << res[i] << "  ";
	if (f)
		cout << endl << setprecision(5) << "\n\tZ_MAX = " << Matrix[n - 1][m - 1] << endl;
	else
		cout << endl << setprecision(5) << "\n\tZ_MIN = " << Matrix[n - 1][m - 1] << endl;

}

void M_method(double**Matrix, double**Matrix1, int n, int m, bool&f1, bool&f2)
{
	int i; bool var = 0, f = 1;
	double*res = inputArray(m);
	cout << "               *******************MINIMUM************************               " << endl;
	res = prepareResult(Matrix1, n, m);
	for (i = 0; i<n - 1; i++)
		add(Matrix1, n, m, n - 1, M, i);
	while (!isMinimumOptimal(Matrix1, n, m) && f){
		f = Simplex_Method(Matrix1, n, m, var);
	}
	if (f){
		cout << endl;
		printSimplexTable(Matrix1, n, m);
		res = prepareResult(Matrix1, n, m);
		if (isArtificialVariableInBasis(res, n, m)){
			cout << "\n No decisions since in basis were artificial variables";
			f1 = 0;
		}
		else
			printResulfLine(Matrix1, res, n, m, 0);
	}
	else {
		cout << endl << "No decisions since all the negative elements in a column" << endl;
		f1 = 0;
	}
	var = 1; f = 1;
	if (debug)
		system("pause");
	cout << endl << "               *******************MAXIMUM************************               " << endl;
	res = prepareResult(Matrix, n, m);
	for (i = 0; i<n - 1; i++)
		add(Matrix, n, m, n - 1, -M, i);
	while (!isMaximumOptimal(Matrix, n, m) && f){
		f = Simplex_Method(Matrix, n, m, var);
	}
	if (f){
		cout << endl;
		printSimplexTable(Matrix, n, m);
		cout << endl;
		res = prepareResult(Matrix, n, m);
		if (isArtificialVariableInBasis(res, n, m)){
			cout << "\n No decisions since in basis were artificial variables";
			f2 = 0;
		}
		else
			printResulfLine(Matrix, res, n, m, 1);
	}
	else {
		cout << endl << "No decisions since all the negative elements in a column" << endl;
		f2 = 0;
	}
	if (debug)
		system("pause");
}


double** deleteArtificialValuses(double**matr, int n, int m)
{
	double**mx = inputMatrix(n, m - n + 1);
	for (int i = 0; i<n; i++)
	for (int j = 0; j<m - n; j++)
		mx[i][j] = matr[i][j];
	for (int i = 0; i<n; i++)
		mx[i][m - n] = matr[i][m - 1];
	return mx;
}


double** Branch_and_Bound(double**Matrix, int n, int m, int*y, int k, bool flag, bool myvalue){
	int i, val; double max = INT_MIN, numb; bool f, f1;
	double**qmatr = inputMatrix(n, m - n + 1);
	qmatr = deleteArtificialValuses(Matrix, n, m);
	m = m - n + 1; int l = 0;
	double*res = inputArray(m);
	res = prepareResult(qmatr, n, m);
	f = 1;
	for (i = 0; i<def; i++)
	{
		f1 = 0;
		for (int j = 0; j<k; j++)
		if (i == y[j] - 1) { f1 = 1; break; }
		if (f1)
		if (fractionalPart(res[i])>pow((double)10, -5) && fractionalPart(res[i])>max) {
			val = i;
			max = fractionalPart(res[i]);
		}
	}
	if (flag) numb = floor(res[val]);
	else numb = ceil(res[val]);
	double**a1 = inputMatrix(n + 1, m + n + 1);
	a1 = deployNewMatrix(qmatr, n, m, numb, val, flag, myvalue);
	n += 1; m += n;
	max = INT_MIN;
	if (!myvalue){
		for (i = 0; i<n - 1; i++)
			add(a1, n, m, n - 1, M, i);
		while (!isMinimumOptimal(a1, n, m) && f)
			f = Simplex_Method(a1, n, m, myvalue);
	}
	else {
		for (i = 0; i<n - 1; i++)
			add(a1, n, m, n - 1, -M, i);
		f = 1;
		while (!isMaximumOptimal(a1, n, m) && f)
			f = Simplex_Method(a1, n, m, myvalue);
	}
	res = prepareResult(a1, n, m);
	if (isArtificialVariableInBasis(res, n, m)) f = 0;
	if (!f) return NULL;
	return a1;

}

double**inputT(double**matr, int n, int m)
{
	double**mx = inputMatrix(n, m);
	equateMatrix(matr, mx, n, m);
	return mx;
}

void Obr(BranchAndBound*&t, double**M, int n, int m, bool flag, bool var, int*y, int k)
{
	if (!t)
	{
		t = new BranchAndBound;
		double*res = inputArray(m);
		res = prepareResult(M, n, m);
		if (!isResulfInteger(res, m, y, k) && Branch_and_Bound(M, n, m, y, k, flag, var) != NULL){
			t->Mx = inputT(Branch_and_Bound(M, n, m, y, k, flag, var), n + 1, m + 2);
			t->n = n + 1;
			t->m = m + 2;
			t->f = 2;
			res = prepareResult(t->Mx, t->n, t->m);
			cout << endl << endl << "Decision: ";
			for (int i = 0; i<def; i++)
				cout << setprecision(5) << res[i] << "   ";
			cout << endl;
			cout << setprecision(5) << "Z = " << t->Mx[t->n - 1][t->m - 1];
			cout << endl << endl;
			if (isResulfInteger(res, t->m, y, k)) t->f = 1;
			else if (Branch_and_Bound(t->Mx, t->n, t->m, y, k, flag, var) == NULL&&Branch_and_Bound(t->Mx, t->n, t->m, y, k, !flag, var) == NULL) t->f = 0;
			t->l = t->r = 0;
		}
		else t = 0;
	}
	else {
		Obr(t->l, t->Mx, t->n, t->m, flag, var, y, k);
		Obr(t->r, t->Mx, t->n, t->m, !flag, var, y, k);
	}
}

bool v_storonu(BranchAndBound*&t)
{
	bool fl = 1;
	if (t)
	{
		if (t->l == 0 && t->r == 0)
		if (t->f == 1 || t->f == 0)
			fl = 1;
		else
			fl = 0;
		return fl&&v_storonu(t->l) && v_storonu(t->r);
	}
	else return fl;
}

bool v_storonu2(BranchAndBound*&t)
{
	bool fl = 1;
	if (t)
	{
		if (t->l == 0 && t->r == 0)
		if (t->f == 0)
			fl = 1;
		else
			fl = 0;
		return fl&&v_storonu2(t->l) && v_storonu2(t->r);
	}
	else return fl;
}

void findMinimumFromBB(BranchAndBound*t, BranchAndBound*&pmin){
	if (t) {
		if (debug) cout << endl << "Element : (" << t->Mx[t->n - 1][t->m - 1] << ") resolution  -->  ";
		if (t->f == 0){
			if (debug) cout << "was zonded" << endl;
		}
		if (t->f == 2)
		if (debug) cout << "incorrectly" << endl;
		if (t->f == 1 && !pmin) {
			pmin = t;
			if (debug) cout << "first $pmin candidat" << endl;
		}
		if (t->f == 1 && t->Mx[t->n - 1][t->m - 1]<pmin->Mx[pmin->n - 1][pmin->m - 1]) {
			pmin = t;
			if (debug) cout << "improves $pmin --> update" << endl;
		}
		findMinimumFromBB(t->l, pmin);
		findMinimumFromBB(t->r, pmin);
	}
}

void findMaximumFromBB(BranchAndBound*t, BranchAndBound*&pmax){
	if (t) {
		if (debug) cout << endl << "Element : (" << t->Mx[t->n - 1][t->m - 1] << ") resolution  -->  ";
		if (t->f == 0){
			if (debug) cout << "was zonded" << endl;
		}
		if (t->f == 2)
		if (debug) cout << "incorrectly" << endl;
		if (t->f == 1 && !pmax) {
			pmax = t;
			if (debug) cout << "first $pmin candidat" << endl;
		}
		if (t->f == 1 && t->Mx[t->n - 1][t->m - 1]>pmax->Mx[pmax->n - 1][pmax->m - 1]) {
			pmax = t;
			if (debug) cout << "improves $pmin --> update" << endl;
		}
		findMaximumFromBB(t->l, pmax);
		findMaximumFromBB(t->r, pmax);
	}
}



void Pechatuska(BranchAndBound*&t, int k){
	/*
	if (t){
	Pechatuska(t->r, k + 3);
	for (int i = 0; i<k; i++)
	cout << " ";
	cout << setprecision(4) << t->Mx[t->n - 1][t->m - 1] << "," << t->f;
	cout << endl << endl;
	Pechatuska(t->l, k + 3);
	}
	*/
	if (t){
		cout << endl << "Domain: " << t->Mx[t->n - 1][t->m - 1] << "(" << t->f << ");" << endl;
		if (t->r)
			cout << "Child R: " << t->r->Mx[t->r->n - 1][t->r->m - 1] << "(" << t->r->f << ");" << endl;
		else
			cout << "Child R: NO" << endl;
		if (t->l)
			cout << "Child L: " << t->l->Mx[t->l->n - 1][t->l->m - 1] << "(" << t->l->f << ");" << endl;
		else
			cout << "Child L: NO" << endl;
		Pechatuska(t->r, k);
		Pechatuska(t->l, k);
	}
}

void show(BranchAndBound*t)
{
	if (t)
	{
		double*res = inputArray(t->m);
		if (t->f == 1)
		{
			res = prepareResult(t->Mx, t->n, t->m);
			for (int i = 0; i < def; i++)
				cout << setprecision(5) << res[i] << "  ";
			cout << "  z=" << t->Mx[t->n - 1][t->m - 1] << endl;
		}
		show(t->l);
		show(t->r);

	}
}

void show1(BranchAndBound*t)
{
	if (t)
	{
		double*res = inputArray(t->m);
		res = prepareResult(t->Mx, t->n, t->m);
		for (int i = 0; i < def; i++)
			cout << setprecision(5) << res[i] << "  ";
		cout << "  z=" << t->Mx[t->n - 1][t->m - 1] << endl;
	}
	show(t->l);
	show(t->r);
}

bool zondirovanie(BranchAndBound*&t, double & optimum, bool minOrMax)
{
	bool f1 = true;
	if (minOrMax&&t->Mx[t->n - 1][t->m - 1] > optimum)
	{
		optimum = t->Mx[t->n - 1][t->m - 1];
		f1 = false;
	}
	if (!minOrMax&&t->Mx[t->n - 1][t->m - 1] < optimum)
	{
		optimum = t->Mx[t->n - 1][t->m - 1];
		f1 = false;
	}
	return f1;
}

void findMinimum(double ** Matrix1, int n, int m, int * y, int numb, BranchAndBound*BranchAndBoundMin, BranchAndBound*p1, double optimum){
	BranchAndBoundMin = new BranchAndBound;
	BranchAndBoundMin->Mx = inputT(Matrix1, n, m);
	BranchAndBoundMin->n = n;
	BranchAndBoundMin->m = m;
	BranchAndBoundMin->f = 2;
	BranchAndBoundMin->l = BranchAndBoundMin->r = 0;
	while (!v_storonu(BranchAndBoundMin)){
		if (zondirovanie(BranchAndBoundMin, optimum, 0)){
			Obr(BranchAndBoundMin, BranchAndBoundMin->Mx, BranchAndBoundMin->n, BranchAndBoundMin->m, 1, 0, y, numb);
		}
	}

	if (debug) {
		cout << endl << "               ******************PECHATUSKA**********************               " << endl;
		Pechatuska(BranchAndBoundMin, 0);
	}
	if (v_storonu2(BranchAndBoundMin)) {
		show1(BranchAndBoundMin);
		cout << endl << "No solutions because all elements in the tree negative" << endl;
		return;
	}
	p1 = 0;
	if (debug) cout << "               *******************MINIMUM-DEBUG-LOG**************               " << endl;
	findMinimumFromBB(BranchAndBoundMin, p1);
	cout << endl << endl << "Decision: ";
	double*res = inputArray(p1->m);
	res = prepareResult(p1->Mx, p1->n, p1->m);
	for (int i = 0; i<def; i++)
		cout << setprecision(5) << res[i] << "   ";
	cout << endl;
	cout << setprecision(5) << "Z = " << p1->Mx[p1->n - 1][p1->m - 1];
	cout << endl << endl;
	show(BranchAndBoundMin);
	delete BranchAndBoundMin;
	delete p1;
}

void findMaximum(double ** Matrix, int n, int m, int * y, int numb, BranchAndBound*BranchAndBoundMax, BranchAndBound*p2, double optimum){
	BranchAndBoundMax = new BranchAndBound;
	BranchAndBoundMax->Mx = inputT(Matrix, n, m);
	BranchAndBoundMax->n = n;
	BranchAndBoundMax->m = m;
	BranchAndBoundMax->f = 2;
	BranchAndBoundMax->l = BranchAndBoundMax->r = 0;
	int k = 0;
	while (!v_storonu(BranchAndBoundMax)){
		if (zondirovanie(BranchAndBoundMax, optimum, 1)){
			Obr(BranchAndBoundMax, BranchAndBoundMax->Mx, BranchAndBoundMax->n, BranchAndBoundMax->m, 1, 1, y, numb);
		}
	}
	if (debug) {
		cout << endl << "               ******************PECHATUSKA**********************               " << endl;
		Pechatuska(BranchAndBoundMax, 0);
	}
	if (v_storonu2(BranchAndBoundMax)) {
		show1(BranchAndBoundMax);
		cout << endl << "No solutions because all elements in the tree negative" << endl;
		return;
	}
	p2 = 0;
	if (debug) cout << "               *******************MAXIMUM-DEBUG-LOG**************               " << endl;
	findMaximumFromBB(BranchAndBoundMax, p2);
	cout << endl << endl << "Decision: ";
	double*res = inputArray(p2->m);
	res = prepareResult(p2->Mx, p2->n, p2->m);
	for (int i = 0; i<def; i++)
		cout << setprecision(5) << res[i] << "   ";
	cout << endl;
	cout << setprecision(5) << "Z = " << p2->Mx[p2->n - 1][p2->m - 1];
	cout << endl << endl;
	show(BranchAndBoundMax);
	delete BranchAndBoundMax;
	delete p2;
}

int main(){
	const int n = 4, m = 9;
	const int N = 5;
	double matr_simply[n][m] = {
		{ 2, -5, 7, 0, 0, 1, 0, 0, 8 },
		{ 1, 0, -8, 7, 2, 0, 1, 0, 1 },
		{ -9, 0, 7, 6, -3, 0, 0, 1, 9 },
		{ -3, 1, -4, 6, -4, M, M, M, 0 },
	};


	setlocale(LC_ALL, ".1251");
	double**Matrix, **Matrix1, *res, *res2;
	int numb; bool f1 = 1, f2 = 1, flag = 1;
	Matrix = inputMatrix(n, m);
	Matrix1 = inputMatrix(n, m);

	for (int i = 0; i<n; i++)
	for (int j = 0; j<m; j++)
		Matrix[i][j] = matr_simply[i][j];


	for (int k = 0; k<n; k++)
	for (int j = 0; j<m; j++)
	if (Matrix[k][j] == M) Matrix1[k][j] = -Matrix[k][j];
	else Matrix1[k][j] = Matrix[k][j];
	cout << endl;
	cout << "**********************************M-METHOD**************************************" << endl;
	M_method(Matrix, Matrix1, n, m, f1, f2);
	res = prepareResult(Matrix, n, m);
	res2 = prepareResult(Matrix1, n, m);
	BranchAndBound *BranchAndBoundMin=0, *BranchAndBoundMax=0, *p1, *p2;
	p1 = 0; p2 = 0;
	int*y;
	do {
		if (debug) 
			system("pause");
		cout << endl << "********************************BRANCH&&BOUND***********************************" << endl;
		cout << "Parameter count:  ";
		cin >> numb;
		y = new int[numb];
		for (int k = 0; k<numb; k++){
			cout << k + 1 << " param index:  ";
			cin >> y[k];
		}

		cout << "               *******************MINIMUM************************               " << endl;
		if (f1&&isResulfInteger(res2, m, y, numb))
			cout << endl << "These variables are integers" << endl;
		else if (f1)
		{
			double optimum = INT_MAX;
			findMinimum(Matrix1, n, m, y, numb, BranchAndBoundMin, p1, optimum);
		}
		if (debug)
			system("pause");
		cout << endl << "               *******************MAXIMUM************************               " << endl;
		if (f2&&isResulfInteger(res, m, y, numb))
			cout << endl << "These variables are integers" << endl;
		else if (f2)
		{
			double optimum = INT_MIN;
			findMaximum(Matrix, n, m, y, numb, BranchAndBoundMax, p2, optimum);
		}
		delete y;
	} while (_getch() != 'q');
	_getch();
	return 0;
}



/* MY
const int n = 4, m = 12;
const int N = 5;
double matr_simply[n][m] = {
{ 1, 1, 8, 6, 6, 1, 0, 0,1,0,0, 40 },
{ -3, 7, 2, 1, -8, 0, 1, 0,0,1,0, 10 },
{ 4, 8, -4, 2, 2, 0, 0, 1,0,0,1, 10 },
{ -4, 0, -4, 5, -5,0,0,0,M, M, M, 0 }
};
*/