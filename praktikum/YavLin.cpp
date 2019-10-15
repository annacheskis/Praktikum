#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

static double X(int i, double h, int num_t) {
	return -1 + (i - num_t) * h;
}

static double F(double v) 
{
	return v / 2;
}

double an_solve(double x) {
	if (x <= 0.25)
		return 1. ;
	if (0.25 < x && x <= 0.5)
		return -4 * x + 2.;
	if (x > 0.5)
		return 0;
	return 0;
}

int max_num_x_ (double tau, double h, int k) {
	double tau1 = tau / k;
	double h1 = h / k;
	int num_t = (int)(1 / tau1);
	int num_x = (int)(2 / h1) + 1 + 2 * num_t;
	return num_x;
}

int main1() {
	std::ofstream out;
	out.open("out_y_l.txt");
	double h, tau;
	int num_t, num_x;
	double x;
	double* cur, * next, * tmp;
	double* cur1;
    //double* v1;
	int i, j;
	double max_delta, sum_delta, abs_max_delta, abs_sum_delta, delta;
	int k = 16;

	tau = 0.1;
	h = 0.1;
    num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;

	cur = (double*)malloc(max_num_x_(tau, h, k) * sizeof(double));
	next = (double*)malloc(max_num_x_(tau, h, k) * sizeof(double));
	cur1 = (double*)malloc(max_num_x_(tau, h, k) * sizeof(double));

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
			next[i] = 0.5 * (cur[i] + cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i])) - (tau / h)*(F(cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i]))) - F(cur[i-1] - (tau / h) * (F(cur[i]) - F(cur[i-1])))));
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
		delta = fabs(an_solve(X(i, h, num_t)) - cur[i]);
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
	
	double tau1 = tau / k; 
	double h1 = h / k;
	num_t = (int)(1 / tau1);
	num_x = (int)(2 / h1) + 1 + 2 * num_t;
	

	for (i = 0; i < num_x; i++) {
		x = X(i, h1, num_t);
		if (x <= -0.25)
			cur1[i] = 1;
		if (-0.25 < x && x <= 0)
			cur1[i] = -4 * x;
		if (x > 0)
			cur1[i] = 0;
		next[i] = cur1[i];
	}

	for (j = 1; j < num_t; j++)
	{
		for (i = 1; i < num_x - 1; i++) {
			//v1[i] = cur[i] - (tau / h) * (F(cur[i + 1]) - F(cur[i]));
			next[i] = 0.5 * (cur1[i] + cur1[i] - (tau1 / h1) * (F(cur1[i + 1]) - F(cur1[i])) - (tau1 / h1)*(F(cur1[i] - (tau1 / h1) * (F(cur1[i + 1]) - F(cur1[i]))) - F(cur1[i - 1] - (tau1 / h1) * (F(cur1[i]) - F(cur1[i - 1])))));
		}
		tmp = cur1;
		cur1 = next;
		next = tmp;
	}

	max_delta = 0;
	abs_max_delta = 0;
	sum_delta = 0;
	abs_sum_delta = 0;

	for (i = num_t; i <= num_x - 1 - num_t; i++) {
		delta = fabs(an_solve(X(i, h1, num_t)) - cur1[i]);
		out << X(i, h1, num_t) << ";" << cur1[i] << "                      " << delta << "\n";

	}
	
	num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;

	//std::cout << X(num_t+1, h, num_t) << " " << X((num_t+1)*k, h1, (int)1 / tau1) << "\n";

	for (i = num_t; i <= num_x - num_t - 1; i++)
	{
		delta = fabs(cur[i] - cur1[k*i]);
		if (delta > max_delta) {
			max_delta = delta;
			std::cout << X(i, h, num_t) << " " << X(k*i, h / k, (int)(1 / tau1)) << " " <<cur[i]<<" "<<cur1[k*i]<<" "<< fabs(cur[i] - cur1[k*i]) << "\n";
		}
		
		if (fabs(cur1[k*i]) > abs_max_delta)
			abs_max_delta = fabs(cur1[k*i]);

		sum_delta = sum_delta + delta;
		abs_sum_delta = abs_sum_delta + fabs(cur1[k*i]);
	}
	abs_max_delta = max_delta / abs_max_delta;
	abs_sum_delta = sum_delta / abs_sum_delta;
	sum_delta = sum_delta * h1;
	
	std::cout << max_delta << " " << sum_delta << " " << abs_max_delta << " " << abs_sum_delta << "\n";
	out << max_delta << " & " << sum_delta << " & " << abs_max_delta << " & " << abs_sum_delta << "\n";

	out.close();
	free(next);
	free(cur);
	free(cur1);
	return 0;
}
