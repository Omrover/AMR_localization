/*
 * Pruebas en host de amr_odom (omni + EKF). gcc, sin RTOS:
 *   bash tests/run_host_tests.sh
 */
#include "ekf.h"
#include "omni.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures;

#define CHECK(cond)                                                               \
	do {                                                                      \
		if (!(cond)) {                                                    \
			printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);    \
			failures++;                                               \
		}                                                                 \
	} while (0)
#define NEAR(a, b, tol) (fabs((double)(a) - (double)(b)) <= (tol))

static const float ANGLES[4] = {0.785398163f, 2.35619449f, 3.92699082f, 5.49778714f};

/* Gaussiana reproducible (Box-Muller con LCG propio). */
static unsigned long rng_state = 12345UL;
static double urand(void)
{
	rng_state = rng_state * 1103515245UL + 12345UL;
	return ((double)((rng_state >> 8) & 0xFFFFFF) + 0.5) / 16777216.0;
}
static float gauss(float std)
{
	double u1 = urand(), u2 = urand();
	return (float)(std * sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2));
}

static void test_omni(void)
{
	amr_omni_t k;
	float w[4], vx, vy, wz;
	const float s = 0.70710678f / 0.05f;

	CHECK(amr_omni_init(&k, 0.05f, 0.25f, ANGLES));
	/* Mismos valores que test_units.py (bridge Python). */
	amr_omni_inverse(&k, 1.0f, 0.0f, 0.0f, w);
	CHECK(NEAR(w[0], -s, 1e-4) && NEAR(w[1], -s, 1e-4) && NEAR(w[2], s, 1e-4) && NEAR(w[3], s, 1e-4));
	amr_omni_inverse(&k, 0.0f, 0.0f, 1.0f, w);
	for (int i = 0; i < 4; i++) {
		CHECK(NEAR(w[i], 5.0, 1e-5));
	}
	for (int n = 0; n < 200; n++) {
		float v[3] = {(float)urand() * 2 - 1, (float)urand() * 2 - 1, (float)urand() * 4 - 2};
		amr_omni_inverse(&k, v[0], v[1], v[2], w);
		amr_omni_forward(&k, w, &vx, &vy, &wz);
		CHECK(NEAR(vx, v[0], 1e-5) && NEAR(vy, v[1], 1e-5) && NEAR(wz, v[2], 1e-5));
		CHECK(amr_omni_residual(&k, w) < 1e-5f);
	}
	/* Una rueda patinando produce residuo. */
	amr_omni_inverse(&k, 0.5f, 0.0f, 0.0f, w);
	w[0] *= 0.5f;
	CHECK(amr_omni_residual(&k, w) > 0.02f);
	/* Configuración inválida. */
	CHECK(!amr_omni_init(&k, 0.0f, 0.25f, ANGLES));
	{
		const float bad[4] = {0, 0, 0, 0};
		CHECK(!amr_omni_init(&k, 0.05f, 0.25f, bad));
	}
}

static void test_utils(void)
{
	CHECK(NEAR(amr_wrap_pi(3.5f), 3.5 - 2 * M_PI, 1e-5));
	CHECK(NEAR(amr_wrap_pi(-3.5f), -3.5 + 2 * M_PI, 1e-5));
	/* yaw de 90 grados alrededor de z */
	CHECK(NEAR(amr_quat_yaw(0.70710678f, 0, 0, 0.70710678f), M_PI / 2, 1e-5));
}

typedef struct {
	double x, y, yaw;
} pose_t;

/*
 * Trayectoria: avanza en curva con traslación lateral, 60 s a 50 Hz.
 * Ruedas con error de escala del 3% en una rueda (desgaste) y ruido;
 * gyro con ruido; yaw de la BNO080 con ruido y una deriva lenta.
 */
static void simulate(bool use_imu, pose_t *final_truth, pose_t *final_est, uint32_t *rejected)
{
	amr_omni_t k;
	amr_ekf_t f;
	amr_ekf_cfg_t cfg;
	const float dt = 0.02f;
	pose_t t = {0, 0, 0};
	const float imu_offset = 1.0f; /* el GRV arranca con un yaw arbitrario */

	amr_omni_init(&k, 0.05f, 0.25f, ANGLES);
	amr_ekf_default_cfg(&cfg);
	amr_ekf_init(&f, &cfg);
	rng_state = 4242UL;
	for (int n = 0; n < 3000; n++) {
		float tt = n * dt;
		float vx = 0.4f, vy = 0.1f * sinf(0.2f * tt), wz = 0.25f * sinf(0.1f * tt);
		float w[4], mvx, mvy, mwz;

		/* verdad (integración fina) */
		for (int sub = 0; sub < 10; sub++) {
			double h = dt / 10.0;
			t.x += (vx * cos(t.yaw) - vy * sin(t.yaw)) * h;
			t.y += (vx * sin(t.yaw) + vy * cos(t.yaw)) * h;
			t.yaw += wz * h;
		}
		amr_omni_inverse(&k, vx, vy, wz, w);
		w[2] *= 1.03f; /* rueda rr con 3% de error de escala */
		for (int i = 0; i < 4; i++) {
			w[i] += gauss(0.05f);
		}
		amr_omni_forward(&k, w, &mvx, &mvy, &mwz);

		amr_ekf_predict(&f, dt);
		amr_ekf_update_wheels(&f, mvx, mvy, mwz);
		if (use_imu) {
			float yaw_imu = amr_wrap_pi((float)t.yaw - imu_offset + 0.0005f * tt + gauss(0.005f));
			amr_ekf_update_gyro(&f, wz + gauss(0.01f));
			amr_ekf_update_imu_yaw(&f, yaw_imu);
		}
	}
	*final_truth = t;
	final_est->x = f.s[AMR_EKF_X];
	final_est->y = f.s[AMR_EKF_Y];
	final_est->yaw = f.s[AMR_EKF_YAW];
	*rejected = f.rejected;
}

static double pos_err(pose_t a, pose_t b)
{
	return hypot(a.x - b.x, a.y - b.y);
}

static void test_ekf_fusion(void)
{
	pose_t truth, wheels_only, fused;
	uint32_t rej1, rej2;
	double e_w, e_f, yaw_w, yaw_f;

	simulate(false, &truth, &wheels_only, &rej1);
	simulate(true, &truth, &fused, &rej2);
	e_w = pos_err(truth, wheels_only);
	e_f = pos_err(truth, fused);
	yaw_w = fabs(amr_wrap_pi((float)(truth.yaw - wheels_only.yaw)));
	yaw_f = fabs(amr_wrap_pi((float)(truth.yaw - fused.yaw)));
	printf("recorrido ~24 m: error pos ruedas=%.3f m, fusion=%.3f m; yaw ruedas=%.4f, fusion=%.4f rad\n",
	       e_w, e_f, yaw_w, yaw_f);
	CHECK(e_f < e_w);
	CHECK(yaw_f < yaw_w);
	CHECK(yaw_f < 0.05); /* deriva del GRV 0.03 rad en 60 s + ruido */
	CHECK(e_f < 0.5);
	CHECK(rej2 < 30);
}

static void test_reset_and_offset(void)
{
	amr_ekf_t f;
	amr_ekf_cfg_t cfg;

	amr_ekf_default_cfg(&cfg);
	amr_ekf_init(&f, &cfg);
	amr_ekf_reset(&f, 1.0f, 2.0f, 0.5f);
	CHECK(amr_ekf_update_imu_yaw(&f, -2.0f)); /* fija el offset sin mover el estado */
	CHECK(NEAR(f.s[AMR_EKF_YAW], 0.5, 1e-6));
	/* La IMU gira 0.1 rad: el yaw del filtro sigue (cruce de +-pi incluido). */
	for (int i = 0; i < 50; i++) {
		amr_ekf_predict(&f, 0.02f);
		amr_ekf_update_imu_yaw(&f, -2.0f + 0.1f);
	}
	CHECK(NEAR(f.s[AMR_EKF_YAW], 0.6, 0.01));
	CHECK(NEAR(f.s[AMR_EKF_X], 1.0, 1e-3) && NEAR(f.s[AMR_EKF_Y], 2.0, 1e-3));

	amr_ekf_reset(&f, 0, 0, 3.1f);
	amr_ekf_update_imu_yaw(&f, 0.0f);
	for (int i = 0; i < 50; i++) {
		amr_ekf_predict(&f, 0.02f);
		amr_ekf_update_imu_yaw(&f, 0.1f); /* 3.1 + 0.1 cruza pi */
	}
	CHECK(NEAR(amr_wrap_pi(f.s[AMR_EKF_YAW] - 3.2f), 0.0, 0.01));
}

static void test_gate_rejects_slip(void)
{
	amr_ekf_t f;
	amr_ekf_cfg_t cfg;
	uint32_t before;

	amr_ekf_default_cfg(&cfg);
	amr_ekf_init(&f, &cfg);
	/* Robot quieto confirmado por gyro y ruedas durante 1 s. */
	for (int i = 0; i < 50; i++) {
		amr_ekf_predict(&f, 0.02f);
		amr_ekf_update_wheels(&f, 0, 0, 0);
		amr_ekf_update_gyro(&f, 0);
	}
	before = f.rejected;
	/* Una rueda patina violentamente: wz de ruedas 2 rad/s, gyro dice 0. */
	amr_ekf_predict(&f, 0.02f);
	amr_ekf_update_gyro(&f, 0.0f);
	amr_ekf_update_wheels(&f, 0.0f, 0.0f, 2.0f);
	CHECK(f.rejected > before);
	CHECK(fabsf(f.s[AMR_EKF_WZ]) < 0.1f);
	CHECK(!amr_ekf_update_gyro(&f, NAN));
}

static void test_variances_positive(void)
{
	amr_ekf_t f;
	amr_ekf_cfg_t cfg;
	float var[6];

	amr_ekf_default_cfg(&cfg);
	amr_ekf_init(&f, &cfg);
	for (int i = 0; i < 500; i++) {
		amr_ekf_predict(&f, 0.02f);
		amr_ekf_update_wheels(&f, 0.3f, 0.0f, 0.1f);
	}
	amr_ekf_variances(&f, var);
	for (int i = 0; i < 6; i++) {
		CHECK(var[i] > 0.0f && isfinite(var[i]));
	}
	CHECK(var[AMR_EKF_X] > 1e-6f); /* la incertidumbre de pose crece sin referencia absoluta */
}

int main(void)
{
	test_omni();
	test_utils();
	test_ekf_fusion();
	test_reset_and_offset();
	test_gate_rejects_slip();
	test_variances_positive();
	if (failures == 0) {
		printf("amr_odom: PASSED (0 failures)\n");
		return 0;
	}
	printf("amr_odom: FAILED (%d failures)\n", failures);
	return 1;
}
