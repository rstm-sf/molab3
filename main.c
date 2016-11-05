#include "main.h"

point3d_t methodExternalPenalty(point3d_t x0, double r0, double C, double epsilon);
point3d_t methodInternalPenalty(point3d_t x0, double r0, double C, double epsilon);

int8_t main() {
	const point3d_t x0 = { 1, 1, 1 };
	const double    r0 = 0.1;
	const double    C = 2.0;
	const double epsilon = 1e-6;

	point3d_t xmin;
	const uint8_t test = 1; // 1 - External Penalty metod
							// 2 - Internal Penalty metod
	switch (test) {
	case 1:
		xmin = methodExternalPenalty(x0, r0, C, epsilon);
		break;
	case 2:
		xmin = methodInternalPenalty(x0, r0, C, epsilon);
		break;
	default:
		printf("Error! value test\n");
		system("pause");
		return -1;
	}

	printf("\nf((%.2f, %.2f, %.2f)) = %.2f\n", xmin.x[0], xmin.x[1], xmin.x[2], f(xmin));

	system("pause");
	return 0;
}

point3d_t methodExternalPenalty(const point3d_t x0, const double r0, const double C, const double epsilon) {
	printf("Start External Penalty metod...\n");
	printf("x0 = (%.2f, %.2f, %.2f)\n", x0.x[0], x0.x[1], x0.x[2]);
	printf("r0 = %f; C = %f; epsilon = %.2e\n", r0, C, epsilon);
	const uint32_t maxIter = 10000;
	double rk = r0;
	point3d_t xk = x0;
	double pk = 0.0;
	uint32_t k = 0;

	for (k; k < maxIter; ++k) {
		xk = nonlinearConjugateGradientMethod(xk, epsilon, maxIter, rk, p_ext, grad_f_ext);
		pk = p_ext(xk, rk);
		if (pk <= epsilon) {
			break;
		}
		rk *= C;
	}

	printf("Iters = %" PRIu32 "\nP(x, rk) = %.2e\n", k, pk);
	printf("End External Penalty metod.\n");
	return xk;
}

point3d_t methodInternalPenalty(const point3d_t x0, const double r0, const double C, const double epsilon) {
	printf("Start Internal Penalty metod...\n");
	printf("x0 = (%.2f, %.2f, %.2f)\n", x0.x[0], x0.x[1], x0.x[2]);
	printf("r0 = %f; C = %f; epsilon = %.2e\n", r0, C, epsilon);
	const uint32_t maxIter = 100;
	double rk = r0;
	point3d_t xk = x0;
	double pk = 0.0;
	uint32_t k = 0;

	for (k; k <maxIter; ++k) {
		xk = nonlinearConjugateGradientMethod(xk, epsilon, maxIter, rk, p_int, grad_f_int);
		pk = p_int(xk, rk);
		if (fabs(pk) <= epsilon) {
			break;
		}
		rk /= C;
	}

	printf("Iters = %" PRIu32 "\nP(x, rk) = %.2e\n", k, pk);
	printf("End Internal Penalty metod.\n");
	return xk;
}
