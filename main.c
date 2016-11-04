#include "main.h"

point3d_t methodExternalPenalty(point3d_t x0, double r0, double C, double epsilon);

int8_t main() {
	const point3d_t x0   = { 100, 550, 550 };
	const double    r0   = 0.1;
	const double    C    = 2.0;
	const double epsilon = 1e-2;

	const point3d_t x = methodExternalPenalty(x0, r0, C, epsilon);
	printf("\nx1   = %.2f\nx2   = %.2f\nx3   = %.2f\nf(x) = %.2f\n", x.x[0], x.x[1], x.x[2], f(x));

	system("pause");
	return 0;
}

point3d_t methodExternalPenalty(const point3d_t x0, const double r0, const double C, const double epsilon) {
	printf("Start Penalty metod...\n");
	printf("x0 = (%.2f, %.2f, %.2f)\n", x0.x[0], x0.x[1], x0.x[2]);
	printf("r0 = %f; C = %f; epsilon = %.2e\n", r0, C, epsilon);
	const uint32_t maxIter = 100;
	double rk    = r0;
	point3d_t xk = x0;
	double pk    = 0.0;
	uint32_t k   = 0;

	for (k; k < maxIter; ++k) {
		xk = methodNewtonRaphson(xk, epsilon, maxIter, rk, p_ext, hessian_ext, grad_f_ext, findT_ext);
		pk = p_ext(xk, rk);
		if (pk <= epsilon) {
			break;
		}
		rk *= C;
	}

	printf("Iters = %" PRIu32 "\nP(x, rk) = %.2e\n", k, pk);
	printf("End Penalty metod.\n");
	return xk;
}
