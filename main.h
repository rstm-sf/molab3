#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

typedef struct point3d {
	double x[3];
} point3d_t;

typedef struct mat3d {
	double a[9];
} mat3d_t;

double f(point3d_t x);
double g1(point3d_t x);
double g2(point3d_t x);
double p(point3d_t x, double rk);
double helpF(point3d_t x, double rk);
mat3d_t hessian(point3d_t x, double rk);
double dotProduct3d(point3d_t x, point3d_t y);
point3d_t add_vec3d(point3d_t x, point3d_t y);
point3d_t scalarmul_vec3d(double alpha, point3d_t x);
double minor2d(double m1, double m2, double m3, double m4);
double det_mat3d(mat3d_t mat);
mat3d_t inversed_mat3d(mat3d_t mat);
bool isPositiveDefMat3d(mat3d_t mat);
point3d_t grad_f(point3d_t x, double rk);
point3d_t methodNewtonRaphson(const point3d_t x0, const double epsilon, const uint32_t maxIter, double rk);



inline double f(const point3d_t x) {
	double tmp = -13.0;
	for (uint8_t i = 0; i < 3; ++i)
		tmp += (x.x[i] - 1.0) * (x.x[i] - 1.0);
	return tmp;
}

inline double g1(const point3d_t x) {
	const double tmp = x.x[0] * x.x[0] + x.x[1] * x.x[1] + x.x[2] * x.x[2] - 4.0;
	return  tmp <= 0.0 ? 0.0 : tmp;
}

inline double g2(const point3d_t x) {
	const double tmp = -x.x[2];
	return  tmp <= 0.0 ? 0.0 : tmp;
}

inline double p(const point3d_t x, const double rk) {
	const double g1_x = g1(x);
	const double g2_x = g2(x);
	return rk * 0.5 * (g1_x * g1_x + g2_x * g2_x);
}

inline double helpF(const point3d_t x, const double rk) {
	const double g1_x = g1(x);
	const double g2_x = g2(x);
	return f(x) + p(x, rk);
}

inline mat3d_t hessian(const point3d_t x, const double rk) {
	mat3d_t H;
	const double g1_x = g1(x);
	const double g2_x = g2(x);
	if (g1_x > 0.0) {
		const double tmp = g2_x > 0.0 ? -3.0 : -3.5;
		H.a[0] = 2.0 * (x.x[0] + rk * (3.0 * x.x[0] * x.x[0] + x.x[1] * x.x[1] + x.x[2] * x.x[2] + tmp));
		H.a[1] = 4.0 * rk * x.x[0] * x.x[1];
		H.a[2] = 4.0 * rk * x.x[0] * x.x[2];
		H.a[3] = H.a[1];
		H.a[4] = 2.0 * x.x[1] * (1.0 + 2.0 * rk * x.x[1]);
		H.a[5] = 4.0 * rk * x.x[1] * x.x[2];
		H.a[6] = H.a[2];
		H.a[7] = H.a[5];
		H.a[8] = 2.0 * x.x[2] * (1.0 + 2.0 * rk * x.x[2]);
	} else {
		const double tmp = 2.0 * x.x[0];
		H.a[0] = g2_x > 0.0 ? tmp * (1.0 + rk) : tmp;
		H.a[4] = 2.0 * x.x[1];
		H.a[8] = 2.0 * x.x[2];
		H.a[1] = H.a[2] = H.a[3] = H.a[5] = H.a[6] = H.a[7] = 0.0;
	}

	return H;
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

inline double minor2d(const double m1, const double m2, const double m3, const double m4) {
	return m1 * m4 - m3 * m2;
}

inline double det_mat3d(const mat3d_t mat) {
	return mat.a[0] * minor2d(mat.a[4], mat.a[5], mat.a[7], mat.a[8]) +
	       mat.a[1] * minor2d(mat.a[3], mat.a[5], mat.a[6], mat.a[8]) +
	       mat.a[2] * minor2d(mat.a[3], mat.a[4], mat.a[6], mat.a[7]);
}

inline mat3d_t inversed_mat3d(const mat3d_t mat) {
	const double det_mat = det_mat3d(mat);
	mat3d_t inv_mat3d = { 0 };
	if (det_mat != 0.0) {
		const double inv_det_mat = 1.0 / det_mat;
		inv_mat3d.a[0] = mat.a[0] * inv_det_mat * minor2d(mat.a[4], mat.a[5], mat.a[7], mat.a[8]);
		inv_mat3d.a[1] = mat.a[3] * inv_det_mat * minor2d(mat.a[1], mat.a[2], mat.a[7], mat.a[8]);
		inv_mat3d.a[2] = mat.a[6] * inv_det_mat * minor2d(mat.a[1], mat.a[2], mat.a[4], mat.a[5]);
		inv_mat3d.a[3] = mat.a[1] * inv_det_mat * minor2d(mat.a[3], mat.a[5], mat.a[6], mat.a[8]);
		inv_mat3d.a[4] = mat.a[4] * inv_det_mat * minor2d(mat.a[0], mat.a[2], mat.a[6], mat.a[8]);
		inv_mat3d.a[5] = mat.a[7] * inv_det_mat * minor2d(mat.a[0], mat.a[2], mat.a[3], mat.a[5]);
		inv_mat3d.a[6] = mat.a[2] * inv_det_mat * minor2d(mat.a[3], mat.a[4], mat.a[6], mat.a[7]);
		inv_mat3d.a[7] = mat.a[5] * inv_det_mat * minor2d(mat.a[0], mat.a[1], mat.a[6], mat.a[7]);
		inv_mat3d.a[8] = mat.a[8] * inv_det_mat * minor2d(mat.a[0], mat.a[1], mat.a[3], mat.a[4]);
	}
	else {
		printf("Determinant = 0!\n");
	}

	return inv_mat3d;
}

inline bool isPositiveDefMat3d(const mat3d_t mat) {
	// Sylvester's criterion
	const double delta1 = mat.a[0];
	const double delta2 = minor2d(mat.a[0], mat.a[1], mat.a[3], mat.a[4]);
	const double delta3 = det_mat3d(mat);
	if (delta1  > 0.0 || delta2  > 0.0 || delta3  > 0.0) {
		return true;
	} else {
		printf("No Positive-definite matrix!\n");
		return false;
	}
}

inline point3d_t grad_f(const point3d_t x, const double rk) {
	const double g1_x = g1(x);
	const double g2_x = g2(x);
	point3d_t grad;
	const double tmp = g2_x > 0.0 ? 2.0 * x.x[0] : 0.0;

	if (g1_x > 0.0) {
		grad.x[0] = 2.0 * (x.x[0] - 1.0) + rk * (x.x[0] * 2.0 * g1_x + tmp);
		grad.x[1] = 2.0 * (x.x[1] - 1.0) + rk * 2.0 * x.x[1] * g1_x;
		grad.x[2] = 2.0 * (x.x[2] - 1.0) + rk * 2.0 * x.x[2] * g1_x;
	} else {
		grad.x[0] = 2.0 * (x.x[0] - 1.0) + rk * tmp;
		grad.x[1] = 2.0 * (x.x[1] - 1.0);
		grad.x[2] = 2.0 * (x.x[2] - 1.0);
	}

	return grad;
}

point3d_t methodNewtonRaphson(const point3d_t x0, const double epsilon, const uint32_t maxIter, const double rk) {
	const double epsilon2 = epsilon * epsilon;
	point3d_t xk1 = x0;
	bool is_seq = false;
	uint32_t k = 0;
	double convergence;

	for (k; k < maxIter; ++k) {
		const point3d_t grad_fun_xk1 = grad_f(xk1, rk);
		const double    dot_grad_xk1 = dotProduct3d(grad_fun_xk1, grad_fun_xk1);
		if (dot_grad_xk1 < epsilon2) {
			printf("grad f(xk) = %.e < epsilon = %.2e\n", sqrt(dot_grad_xk1), epsilon);
			break;
		}

		const mat3d_t H = hessian(xk1, rk);
		const mat3d_t invH = inversed_mat3d(H);
		const point3d_t d = isPositiveDefMat3d(invH) ? mat_vec3d(-1.0, H, grad_fun_xk1)
		                                             : scalarmul_vec3d(-1.0, grad_fun_xk1);

		const point3d_t xk = xk1;
		const double t = findT_phi(xk, d);
		const point3d_t xk1_minus_xk = scalarmul_vec3d(t, d);
		xk1 = add_vec3d(xk1, xk1_minus_xk);

		convergence = dotProduct2d(xk1_minus_xk, xk1_minus_xk);
		const double abs_fk1_minus_fk = fabs(fun(xk1) - fun(xk));
		if (convergence < epsilon2 && abs_fk1_minus_fk < epsilon) {
			if (is_seq) {
				break;
			} else {
				is_seq = true;
			}
		} else {
			is_seq = false;
		}
	}

	printf("Iters = %" PRIu32 "\nConvergence = %.2e\n", k, sqrt(convergence));

	return xk1;
}

#endif MAIN_H