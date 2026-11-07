#include "ekf.h"

#include <math.h>
#include <string.h>

#define PI_F 3.14159265358979f

float amr_wrap_pi(float a)
{
	return atan2f(sinf(a), cosf(a));
}

float amr_quat_yaw(float w, float x, float y, float z)
{
	return atan2f(2.0f * (w * z + x * y), 1.0f - 2.0f * (y * y + z * z));
}

void amr_ekf_default_cfg(amr_ekf_cfg_t *cfg)
{
	cfg->accel_noise = 1.0f;
	cfg->alpha_noise = 2.0f;
	cfg->wheel_lin_std = 0.02f;
	cfg->wheel_ang_std = 0.05f;
	cfg->gyro_std = 0.01f;
	cfg->yaw_std = 0.01f;
	cfg->gate_chi2 = 16.0f;
}

void amr_ekf_init(amr_ekf_t *f, const amr_ekf_cfg_t *cfg)
{
	memset(f, 0, sizeof(*f));
	f->cfg = *cfg;
	amr_ekf_reset(f, 0.0f, 0.0f, 0.0f);
}

void amr_ekf_reset(amr_ekf_t *f, float x, float y, float yaw)
{
	memset(f->s, 0, sizeof(f->s));
	memset(f->P, 0, sizeof(f->P));
	f->s[AMR_EKF_X] = x;
	f->s[AMR_EKF_Y] = y;
	f->s[AMR_EKF_YAW] = amr_wrap_pi(yaw);
	f->P[AMR_EKF_X][AMR_EKF_X] = 1e-6f;
	f->P[AMR_EKF_Y][AMR_EKF_Y] = 1e-6f;
	f->P[AMR_EKF_YAW][AMR_EKF_YAW] = 1e-6f;
	f->P[AMR_EKF_VX][AMR_EKF_VX] = 1.0f;
	f->P[AMR_EKF_VY][AMR_EKF_VY] = 1.0f;
	f->P[AMR_EKF_WZ][AMR_EKF_WZ] = 1.0f;
	f->yaw_offset_valid = false;
}

static void symmetrize(amr_ekf_t *f)
{
	for (int i = 0; i < AMR_EKF_N; i++) {
		for (int j = i + 1; j < AMR_EKF_N; j++) {
			float m = 0.5f * (f->P[i][j] + f->P[j][i]);
			f->P[i][j] = m;
			f->P[j][i] = m;
		}
		if (f->P[i][i] < 1e-12f) {
			f->P[i][i] = 1e-12f;
		}
	}
}

void amr_ekf_predict(amr_ekf_t *f, float dt)
{
	float F[AMR_EKF_N][AMR_EKF_N] = {{0}};
	float FP[AMR_EKF_N][AMR_EKF_N];
	float c, s, vx, vy;

	if (!(dt > 0.0f)) {
		return;
	}
	c = cosf(f->s[AMR_EKF_YAW]);
	s = sinf(f->s[AMR_EKF_YAW]);
	vx = f->s[AMR_EKF_VX];
	vy = f->s[AMR_EKF_VY];

	for (int i = 0; i < AMR_EKF_N; i++) {
		F[i][i] = 1.0f;
	}
	F[AMR_EKF_X][AMR_EKF_YAW] = (-vx * s - vy * c) * dt;
	F[AMR_EKF_X][AMR_EKF_VX] = c * dt;
	F[AMR_EKF_X][AMR_EKF_VY] = -s * dt;
	F[AMR_EKF_Y][AMR_EKF_YAW] = (vx * c - vy * s) * dt;
	F[AMR_EKF_Y][AMR_EKF_VX] = s * dt;
	F[AMR_EKF_Y][AMR_EKF_VY] = c * dt;
	F[AMR_EKF_YAW][AMR_EKF_WZ] = dt;

	f->s[AMR_EKF_X] += (vx * c - vy * s) * dt;
	f->s[AMR_EKF_Y] += (vx * s + vy * c) * dt;
	f->s[AMR_EKF_YAW] = amr_wrap_pi(f->s[AMR_EKF_YAW] + f->s[AMR_EKF_WZ] * dt);

	for (int i = 0; i < AMR_EKF_N; i++) {
		for (int j = 0; j < AMR_EKF_N; j++) {
			float acc = 0.0f;
			for (int k = 0; k < AMR_EKF_N; k++) {
				acc += F[i][k] * f->P[k][j];
			}
			FP[i][j] = acc;
		}
	}
	for (int i = 0; i < AMR_EKF_N; i++) {
		for (int j = 0; j < AMR_EKF_N; j++) {
			float acc = 0.0f;
			for (int k = 0; k < AMR_EKF_N; k++) {
				acc += FP[i][k] * F[j][k];
			}
			f->P[i][j] = acc;
		}
	}
	f->P[AMR_EKF_VX][AMR_EKF_VX] += f->cfg.accel_noise * f->cfg.accel_noise * dt;
	f->P[AMR_EKF_VY][AMR_EKF_VY] += f->cfg.accel_noise * f->cfg.accel_noise * dt;
	f->P[AMR_EKF_WZ][AMR_EKF_WZ] += f->cfg.alpha_noise * f->cfg.alpha_noise * dt;
	symmetrize(f);
}

/* Corrección escalar de la componente k con H = e_k. */
static bool update_scalar(amr_ekf_t *f, int k, float z, float std, bool angle)
{
	float innov = z - f->s[k];
	float S, K[AMR_EKF_N], row[AMR_EKF_N];

	if (!isfinite(z)) {
		f->rejected++;
		return false;
	}
	if (angle) {
		innov = amr_wrap_pi(innov);
	}
	S = f->P[k][k] + std * std;
	if (innov * innov > f->cfg.gate_chi2 * S) {
		f->rejected++;
		return false;
	}
	for (int i = 0; i < AMR_EKF_N; i++) {
		K[i] = f->P[i][k] / S;
		row[i] = f->P[k][i];
	}
	for (int i = 0; i < AMR_EKF_N; i++) {
		f->s[i] += K[i] * innov;
		for (int j = 0; j < AMR_EKF_N; j++) {
			f->P[i][j] -= K[i] * row[j];
		}
	}
	f->s[AMR_EKF_YAW] = amr_wrap_pi(f->s[AMR_EKF_YAW]);
	symmetrize(f);
	return true;
}

int amr_ekf_update_wheels(amr_ekf_t *f, float vx, float vy, float wz)
{
	int ok = 0;

	ok += update_scalar(f, AMR_EKF_VX, vx, f->cfg.wheel_lin_std, false);
	ok += update_scalar(f, AMR_EKF_VY, vy, f->cfg.wheel_lin_std, false);
	ok += update_scalar(f, AMR_EKF_WZ, wz, f->cfg.wheel_ang_std, false);
	return ok;
}

bool amr_ekf_update_gyro(amr_ekf_t *f, float wz)
{
	return update_scalar(f, AMR_EKF_WZ, wz, f->cfg.gyro_std, false);
}

bool amr_ekf_update_imu_yaw(amr_ekf_t *f, float yaw_imu)
{
	if (!isfinite(yaw_imu)) {
		f->rejected++;
		return false;
	}
	if (!f->yaw_offset_valid) {
		f->yaw_offset = amr_wrap_pi(f->s[AMR_EKF_YAW] - yaw_imu);
		f->yaw_offset_valid = true;
		return true;
	}
	return update_scalar(f, AMR_EKF_YAW, amr_wrap_pi(yaw_imu + f->yaw_offset), f->cfg.yaw_std,
			     true);
}

void amr_ekf_variances(const amr_ekf_t *f, float out[6])
{
	for (int i = 0; i < AMR_EKF_N; i++) {
		out[i] = f->P[i][i];
	}
}
