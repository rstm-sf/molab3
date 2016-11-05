#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>

typedef struct segment {
	double a;
	double b;
} segment_t;

typedef struct point3d {
	double x[3];
} point3d_t;

typedef struct mat3d {
	double a[9];
} mat3d_t;

inline double f(const point3d_t x) {
	double tmp = -13.0;
	for (uint8_t i = 0; i < 3; ++i) {
		const double tmp2 = x.x[i] - 1.0;
		tmp += tmp2 * tmp2;
	}

	return tmp;
}

inline double g1(const point3d_t x) {
	return x.x[0] * x.x[0] + x.x[1] * x.x[1] + x.x[2] * x.x[2] - 4.0;
}

inline double g1_cut(const point3d_t x) {
	const double tmp = g1(x);
	return  tmp <= 0.0 ? 0.0 : tmp;
}

inline double g2(const point3d_t x) {
	return -x.x[2];
}

inline double g2_cut(const point3d_t x) {
	const double tmp = g2(x);
	return  tmp <= 0.0 ? 0.0 : tmp;
}

inline double p_ext(const point3d_t x, const double rk) {
	const double g1_x = g1_cut(x);
	const double g2_x = g2_cut(x);
	return rk * 0.5 * (g1_x * g1_x + g2_x * g2_x);
}

inline double p_int(const point3d_t x, const double rk) {
	const double g1_x = g1(x);
	const double g2_x = g2(x);
	if (g1_x >= 0 || g2_x >= 0) {
		return -(double)INFINITY;
	}
	return -rk * (log(-g1_x) + log(-g2_x));
}

inline double dotProduct3d(const point3d_t x, const point3d_t y) {
	return x.x[0] * y.x[0] + x.x[1] * y.x[1] + x.x[2] * y.x[2];
}

inline point3d_t add_vec3d(const point3d_t x, const point3d_t y) {
	const point3d_t vec = { x.x[0] + y.x[0], x.x[1] + y.x[1], x.x[2] + y.x[2] };
	return vec;
}

inline point3d_t scalarmul_vec3d(const double alpha, const point3d_t x) {
	const point3d_t vec = { alpha * x.x[0], alpha * x.x[1], alpha * x.x[2] };
	return vec;
}

inline point3d_t mat_vec3d(const double alpha, const mat3d_t mat, const point3d_t x) {
	const point3d_t vec = { alpha * (mat.a[0] * x.x[0] + mat.a[1] * x.x[1] + mat.a[2] * x.x[2]),
	                        alpha * (mat.a[3] * x.x[0] + mat.a[4] * x.x[1] + mat.a[5] * x.x[2]),
	                        alpha * (mat.a[6] * x.x[0] + mat.a[7] * x.x[1] + mat.a[8] * x.x[2]) };
	return vec;
}

inline point3d_t grad_f_ext(const point3d_t x, const double rk) {
	const double g1_x = g1_cut(x);
	const double g2_x = g2_cut(x);
	point3d_t grad;
	const double tmp = g2_x > 0.0 ? 2.0 * x.x[2] : 0.0;

	if (g1_x > 0.0) {
		grad.x[0] = 2.0 * (x.x[0] - 1.0) + rk * 2.0 * x.x[0] * g1_x;
		grad.x[1] = 2.0 * (x.x[1] - 1.0) + rk * 2.0 * x.x[1] * g1_x;
		grad.x[2] = 2.0 * (x.x[2] - 1.0) + rk * (x.x[2] * 2.0 * g1_x + tmp);
	} else {
		grad.x[0] = 2.0 * (x.x[0] - 1.0);
		grad.x[1] = 2.0 * (x.x[1] - 1.0);
		grad.x[2] = 2.0 * (x.x[2] - 1.0) + rk * tmp;
	}

	return grad;
}

inline point3d_t grad_f_int(const point3d_t x, const double rk) {
	const double g1_x = -g1(x);
	const double inv_g1_x = 1.0 / g1_x;
	const point3d_t grad = {
		2.0 * (x.x[0] * (1.0 + rk * inv_g1_x) - 1.0),
		2.0 * (x.x[1] * (1.0 + rk * inv_g1_x) - 1.0),
		2.0 * (x.x[2] - 1.0 + rk * (x.x[2] * inv_g1_x - 0.5 / x.x[2]))
	};

	return grad;
}

inline double phi(const double t, const point3d_t x, const double rk, const point3d_t d, double(*p)(point3d_t, double)) {
	const point3d_t xtd = add_vec3d(x, scalarmul_vec3d(t, d));
	return f(xtd) + p(xtd, rk);
}

double methodGoldenSection(const segment_t seg, const double epsilon, const point3d_t x, const double rk, const point3d_t d, double(*p)(point3d_t, double)) {
	//printf("Start The method of the Golden section...\n");
	double ak = seg.a;
	double bk = seg.b;
	// (3 - sqrt(5))/2 ~ 0.3819660112501051
	double yk,zk, f_yk, f_zk;
	do {
		bk = bk / 2;
		yk = ak + 0.3819660112501051 * (bk - ak);
		zk = ak + bk - yk;
		f_yk = phi(yk, x, rk, d, p);
		f_zk = phi(zk, x, rk, d, p);
	} while (f_yk == -(double)INFINITY || f_zk == -(double)INFINITY);
	double convergence;

	uint16_t k = 0;
	do {
		if (f_yk <= f_zk) {
			bk = zk;
			zk = yk;
			f_zk = f_yk;
			yk = ak + bk - yk;
			f_yk = phi(yk, x, rk, d, p);
		} else {
			ak = yk;
			yk = zk;
			f_yk = f_zk;
			zk = ak + bk - zk;
			f_zk = phi(zk, x, rk, d, p);
		}
		convergence = fabs(ak - bk);
		k++;
	} while (convergence > epsilon);

	//printf("End The method of the Golden section\nIters = %" PRIu16 "\nConvergence = %e\n", k, convergence);

	return (ak + bk) * 0.5;
}

inline double findT(const point3d_t x, const double rk, const point3d_t d, double(*p)(point3d_t, double)) {
	const segment_t seg = { 0, 100 };
	const double epsilon = 1.0e-6;
	const double tmin = methodGoldenSection(seg, epsilon, x, rk, d, p);
	return tmin;
}

point3d_t nonlinearConjugateGradientMethod(const point3d_t x0, const double epsilon, const uint32_t maxIter, const double rk,
	double(*p)(point3d_t, double), point3d_t(*grad_fun)(point3d_t, double)) {
	//printf("Start Nonlinear Conjugate Gradient Method...\n");
	const double epsilon2 = epsilon * epsilon;
	point3d_t xk1 = x0;
	point3d_t d;
	bool is_seq = false;
	uint32_t k = 0;
	double dot_grad_xk, convergence;

	for (k; k < maxIter; ++k) {
		const point3d_t grad_fun_xk1 = grad_fun(xk1, rk);
		const double    dot_grad_xk1 = dotProduct3d(grad_fun_xk1, grad_fun_xk1);
		if (dot_grad_xk1 < epsilon2) {
			//printf("grad f(xk) = %.e < epsilon = %.2e\n", sqrt(dot_grad_xk1), epsilon);
			break;
		}

		if (k == 0) {
			d = scalarmul_vec3d(-1.0, grad_fun_xk1);
		} else {
			const double beta = dot_grad_xk1 / dot_grad_xk;
			d = add_vec3d(scalarmul_vec3d(beta, d), scalarmul_vec3d(-1.0, grad_fun_xk1));
		}

		const point3d_t xk = xk1;
		const double t = findT(xk, rk, d, p);
		const point3d_t xk1_minus_xk = scalarmul_vec3d(t, d);
		xk1 = add_vec3d(xk1, xk1_minus_xk);

		convergence = dotProduct3d(xk1_minus_xk, xk1_minus_xk);
		const double abs_fk1_minus_fk = fabs(p(xk1, rk) - p(xk, rk));
		if (convergence < epsilon2 && abs_fk1_minus_fk < epsilon) {
			if (is_seq) {
				break;
			} else {
				is_seq = true;
			}
		} else {
			is_seq = false;
		}
		dot_grad_xk = dot_grad_xk1;
	}

	//printf("Iters = %" PRIu32 "\nConvergence = %.2e\n", k, sqrt(convergence));
	//printf("End Nonlinear Conjugate Gradient Method.\n");
	return xk1;
}

#endif MAIN_H
