#include "main.h"

point3d_t methodPenalty(point3d_t x0, double r0, double C, double epsilon);

int8_t main() {
	printf("Sleep.\n");
	const point3d_t x0   = { 0, 0, 0 };
	const double    r0   = 0.1;
	const double    C    = 2.0;
	const double epsilon = 1e-2;

	const point3d_t x = methodPenalty(x0, r0, C, epsilon);
	printf("x = %f\ny = %f\nz = %f\n", x.x[0], x.x[1], x.x[2]);

	system("pause");
	return 0;
}

point3d_t methodPenalty(const point3d_t x0, const double r0, const double C, const double epsilon) {
	printf("Penalty metod\n");
	const uint32_t maxIter = 100;
	double rk    = r0;
	point3d_t xk = x0;
	double pk    = 0.0;
	uint32_t k   = 0;

	for (k; k < maxIter; ++k) {
		xk = methodNewtonRaphson(xk, epsilon, maxIter, rk);
		pk = p(xk, rk);
		if (pk <= epsilon) {
			break;
		}
		rk *= C;
	}

	printf("Iters = %" PRIu32 "\nP = %.2e\n", k, pk);
	return xk;
}
