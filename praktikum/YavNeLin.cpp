#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

static double X(int i, double h, int num_t)
{
	return -1 + (i - num_t) * h;
}

static double F(double v)
{
    if(v < 0.00000000000001)
        return 0;
    else
        return v * v / 2;
}

double an_solve_nelin(double x) 
{
	    if (x < 0.75)
			return 1;
		if (fabs(x - 0.75) <0.000000001)
			return 0.5;
		if (x > 0.75)
			return 0;
		return 0;
}

int max_num_x__(double tau, double h, int k) {
	double tau1 = tau / k;
	double h1 = h / k;
	int num_t = (int)(1 / tau1);
	int num_x = (int)(2 / h1) + 1 + 2 * num_t;
	return num_x;
}

void calculate2(double h, double tau, double *cur, double *next, std::ofstream &out)
{
	int num_t, num_x;
	double x;
	double *tmp;
	int i, j;
	double max_delta, sum_delta, abs_max_delta, abs_sum_delta, delta;

	num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;


	for (i = 0; i < num_x; i++) {
		x = X(i, h, num_t);
		if (x <= -0.25)
			cur[i] = 1;
		if (-0.25 < x && x <= 0)
			cur[i] = -4 * x;
		if (x > 0)
			cur[i] = 0;
		next[i] = cur[i];
	}

	for (j = 1; j < num_t; j++)
	{
		for (i = 1; i < num_x - 1; i++) {
			next[i] = 0.5 * (cur[i] + cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i])) - (tau / h) * (F(cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i]))) - F(cur[i - 1] - (tau / h) * (F(cur[i]) - F(cur[i - 1])))));
		}
		tmp = cur;
		cur = next;
		next = tmp;
	}


	max_delta = 0;
	sum_delta = 0;
	abs_max_delta = 0;
	abs_sum_delta = 0;

	for (i = num_t; i <= num_x - num_t - 1; i++)
	{
		delta = fabs(an_solve_nelin(X(i, h, num_t)) - cur[i]);
		out << X(i, h, num_t) << ";" << cur[i] << "                      " << delta << "\n";
		if (fabs(cur[i]) >= abs_max_delta)
			abs_max_delta = fabs(cur[i]);
		if (delta >= max_delta)
			max_delta = delta;
		sum_delta += delta;
		abs_sum_delta += cur[i];
	}

	abs_max_delta = max_delta / abs_max_delta;
	sum_delta = h * sum_delta;
	abs_sum_delta = sum_delta / (h * abs_sum_delta);

	std::cout << max_delta << " " << sum_delta << " " << abs_max_delta << " " << abs_sum_delta << "\n";
	out << max_delta << " & " << sum_delta << " & " << abs_max_delta << " & " << abs_sum_delta << "\n";
}

void compare2(double *cur, double *cur1, int k, std::ofstream &out, double h, double tau) {
	int i;
	double max_delta, sum_delta, abs_max_delta, abs_sum_delta, delta;
	int num_x, num_t;

	num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;

	max_delta = 0;
	sum_delta = 0;
	abs_max_delta = 0;
	abs_sum_delta = 0;
	for (i = num_t; i <= num_x - num_t - 1; i++)
	{
		delta = fabs(cur[i] - cur1[k*i]);
		if (delta > max_delta) {
			max_delta = delta;
		}

		if (fabs(cur1[k*i]) > abs_max_delta)
			abs_max_delta = fabs(cur1[k*i]);

		sum_delta = sum_delta + delta;
		abs_sum_delta = abs_sum_delta + fabs(cur1[k*i]);
	}
	abs_max_delta = max_delta / abs_max_delta;
	abs_sum_delta = sum_delta / abs_sum_delta;
	sum_delta = sum_delta * h / k;

	std::cout << max_delta << " " << sum_delta << " " << abs_max_delta << " " << abs_sum_delta << "\n";
	out << max_delta << " & " << sum_delta << " & " << abs_max_delta << " & " << abs_sum_delta << "\n";

}

int main()
{
	std::ofstream out;
	out.open("out_y_nl.txt");
	double h, tau;
	int k = 4;
	tau = 0.01;
	h = 0.01;
	double* cur, *next, *cur1;

	cur = (double*)malloc(max_num_x__(tau, h, k) * sizeof(double));
	next = (double*)malloc(max_num_x__(tau, h, k) * sizeof(double));
	cur1 = (double*)malloc(max_num_x__(tau, h, k) * sizeof(double));

	calculate2(h, tau, cur, next, out);

	calculate2(h / k, tau / k, cur1, next, out);

	compare2(cur, cur1, k, out, h, tau);

	out.close();
	free(next);
	free(cur);
	free(cur1);
	return 0;
}
