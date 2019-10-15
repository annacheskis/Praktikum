#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

static double X(int i, double h, int num_t) {
	return -1 + (i - num_t) * h;
}


double an_solve_ny(double x) {
	if (x <= 0.25)
		return 1.;
	if (0.25 < x && x <= 0.5)
		return -4 * x + 2.;
	if (x > 0.5)
		return 0;
	return 0;
}

int max_num_x(double tau, double h, int k) {
	double tau1 = tau / k;
	double h1 = h / k;
	int num_t = (int)(1 / tau1);
	int num_x = (int)(2 / h1) + 1 + 2 * num_t;
	return num_x;
}

void calculate(double h, double tau, double *cur, double *next, std::ofstream &out)
{
	double c1;
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

	c1 = tau / (4 * h - tau);

	for (j = 1; j < num_t; j++)
	{
		next[num_x - 1] = 1;
		for (i = num_x - 2; i > 0; i--) {
			next[i] = cur[i] + c1 * (cur[i - 1] - next[i + 1]);
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
		delta = fabs(an_solve_ny(X(i, h, num_t)) - cur[i]);
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

int main() {
	std::ofstream out;
	out.open("out_ny_l.txt");
	double h, tau;
	int k = 2;
	tau = 0.01;
	h = 0.01;
	double* cur, *next, *cur1;
	cur = (double*)malloc(max_num_x(tau, h, k) * sizeof(double));
	next = (double*)malloc(max_num_x(tau, h, k) * sizeof(double));
	cur1 = (double*)malloc(max_num_x(tau, h, k) * sizeof(double));

	calculate(h, tau, cur, next, out);

	
	/*num_t = (int)(1 / tau);
	num_x = (int)(2 / h) + 1 + 2 * num_t;

	//std::cout << X(num_t+1, h, num_t) << " " << X((num_t+1)*k, h1, (int)1 / tau1) << "\n";

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
	sum_delta = sum_delta * h1;

	std::cout << max_delta << " " << sum_delta << " " << abs_max_delta << " " << abs_sum_delta << "\n";
	out << max_delta << " & " << sum_delta << " & " << abs_max_delta << " & " << abs_sum_delta << "\n";

	*/
	out.close();
	free(next);
	free(cur);
	free(cur1);
	return 0;
}
