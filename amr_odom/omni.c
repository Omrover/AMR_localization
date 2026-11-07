#include "omni.h"

#include <math.h>

bool amr_omni_init(amr_omni_t *k, float wheel_radius, float center_distance,
		   const float angles_rad[AMR_OMNI_WHEELS])
{
	double a[3][3] = {{0}};
	double inv[3][3];
	double det;

	if (!(wheel_radius > 0.0f) || !(center_distance > 0.0f)) {
		return false;
	}
	k->r = wheel_radius;
	k->L = center_distance;
	for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
		k->J[i][0] = -sinf(angles_rad[i]);
		k->J[i][1] = cosf(angles_rad[i]);
		k->J[i][2] = center_distance;
	}
	/* A = J^T J (3x3) e inversa por cofactores, en double. */
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
				a[r][c] += (double)k->J[i][r] * (double)k->J[i][c];
			}
		}
	}
	det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
	      a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
	      a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
	if (fabs(det) < 1e-9) {
		return false;
	}
	inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
	for (int r = 0; r < 3; r++) {
		for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
			double s = 0.0;
			for (int c = 0; c < 3; c++) {
				s += inv[r][c] * (double)k->J[i][c];
			}
			k->Jp[r][i] = (float)s;
		}
	}
	return true;
}

void amr_omni_inverse(const amr_omni_t *k, float vx, float vy, float wz,
		      float wheel_rad_s[AMR_OMNI_WHEELS])
{
	for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
		wheel_rad_s[i] = (k->J[i][0] * vx + k->J[i][1] * vy + k->J[i][2] * wz) / k->r;
	}
}

void amr_omni_forward(const amr_omni_t *k, const float wheel_rad_s[AMR_OMNI_WHEELS],
		      float *vx, float *vy, float *wz)
{
	float out[3];

	for (int r = 0; r < 3; r++) {
		float s = 0.0f;
		for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
			s += k->Jp[r][i] * wheel_rad_s[i] * k->r;
		}
		out[r] = s;
	}
	*vx = out[0];
	*vy = out[1];
	*wz = out[2];
}

float amr_omni_residual(const amr_omni_t *k, const float wheel_rad_s[AMR_OMNI_WHEELS])
{
	float vx, vy, wz, pred[AMR_OMNI_WHEELS], acc = 0.0f;

	amr_omni_forward(k, wheel_rad_s, &vx, &vy, &wz);
	amr_omni_inverse(k, vx, vy, wz, pred);
	for (int i = 0; i < AMR_OMNI_WHEELS; i++) {
		float e = (wheel_rad_s[i] - pred[i]) * k->r;
		acc += e * e;
	}
	return sqrtf(acc / AMR_OMNI_WHEELS);
}
