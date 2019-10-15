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

int main2()
{
	std::ofstream out;
	out.open("out_y_nl.txt");
	double h, tau;
	int num_t, num_x;
	double x;
	double* cur, * next, * tmp;
	int i, j;
	double max_delta, sum_delta, abs_max_delta, abs_sum_delta, delta;

	tau = 0.001;
	h = 0.01;
	num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;
	cur = (double*)malloc(num_x * sizeof(double));
	next = (double*)malloc(num_x * sizeof(double));


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
			next[i] = 0.5 * (cur[i] + cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i])) - (tau / h) * (F(cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i]))) - F(cur[i-1] - (tau / h) * (F(cur[i]) - F(cur[i-1])))));
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

	out.close();
	free(next);
	free(cur);
	return 0;
}
