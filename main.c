#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

typedef struct point3d {
	double x[3];
} point3d_t;

double f(point3d_t x);
bool g1(point3d_t x);
bool g2(point3d_t x);

point3d_t methodPenalty(point3d_t x0, double r0, double epsilon);

int8_t main() {
	printf("Sleep.\n");

	system("pause");
	return 0;
}

inline double f(const point3d_t x) {
	double tmp = -13.0;
	for (uint8_t i = 0; i < 3; ++i)
		tmp += (x.x[i] - 1.0) * (x.x[i] - 1.0);
	return tmp;
}

inline bool g1(const point3d_t x) {
	return x.x[0] * x.x[0] + x.x[1] * x.x[1] + x.x[2] * x.x[2] <= 4.0 ? true : false;
}

inline bool g2(const point3d_t x) {
	return x.x[2] >= 0.0 ? true : false;
}