/*
 * EKF de odometría del MCU: fusiona velocidades de rueda (encoders vía RoboClaw)
 * con la BNO080 (Game Rotation Vector y gyro calibrado; la BNO080 ya ejecuta su
 * propio filtro de fusión on-chip).
 *
 * Estado s = [x, y, yaw, vx, vy, wz]: pose en odom y velocidades en base_footprint.
 * Predicción: velocidad constante en el frame del robot.
 * Correcciones escalares secuenciales (R diagonal), así no hay inversiones de matriz:
 *   - ruedas: vx, vy, wz (cinemática directa)
 *   - gyro z: wz
 *   - yaw BNO080: yaw + offset (el GRV es relativo; el offset se fija en el
 *     primer dato tras un reset). Innovación envuelta a [-pi, pi].
 * Cada medición pasa una compuerta de Mahalanobis: si las ruedas patinan y
 * contradicen al gyro, la medición se descarta y se cuenta.
 *
 * Orden recomendado por ciclo: predict -> gyro -> yaw IMU -> ruedas. Así la
 * compuerta compara las ruedas contra la referencia más confiable (el gyro).
 */
#ifndef AMR_EKF_H
#define AMR_EKF_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { AMR_EKF_X, AMR_EKF_Y, AMR_EKF_YAW, AMR_EKF_VX, AMR_EKF_VY, AMR_EKF_WZ, AMR_EKF_N };

typedef struct {
	float accel_noise;      /* densidad de aceleración lineal [m/s^2] */
	float alpha_noise;      /* densidad de aceleración angular [rad/s^2] */
	float wheel_lin_std;    /* ruido de vx, vy de ruedas [m/s] */
	float wheel_ang_std;    /* ruido de wz de ruedas [rad/s] */
	float gyro_std;         /* ruido del gyro z [rad/s] */
	float yaw_std;          /* ruido del yaw de la BNO080 [rad] */
	float gate_chi2;        /* umbral de la compuerta por medición escalar (9 = 3 sigma) */
} amr_ekf_cfg_t;

typedef struct {
	amr_ekf_cfg_t cfg;
	float s[AMR_EKF_N];
	float P[AMR_EKF_N][AMR_EKF_N];
	float yaw_offset;       /* yaw_odom = yaw_imu + offset */
	bool yaw_offset_valid;
	uint32_t rejected;      /* mediciones descartadas por la compuerta */
} amr_ekf_t;

void amr_ekf_default_cfg(amr_ekf_cfg_t *cfg);
void amr_ekf_init(amr_ekf_t *f, const amr_ekf_cfg_t *cfg);
/* Fija la pose (odom_reset). El offset del yaw de la IMU se recalcula con el siguiente dato. */
void amr_ekf_reset(amr_ekf_t *f, float x, float y, float yaw);
void amr_ekf_predict(amr_ekf_t *f, float dt);
/* Devuelve el número de componentes aceptadas (0..3). */
int amr_ekf_update_wheels(amr_ekf_t *f, float vx, float vy, float wz);
bool amr_ekf_update_gyro(amr_ekf_t *f, float wz);
bool amr_ekf_update_imu_yaw(amr_ekf_t *f, float yaw_imu);
/* Varianzas (x, y, yaw, vx, vy, wz) para llc/odom.cov. */
void amr_ekf_variances(const amr_ekf_t *f, float out[6]);

/* Utilidades */
float amr_wrap_pi(float a);
/* Yaw de un cuaternión (w, x, y, z). */
float amr_quat_yaw(float w, float x, float y, float z);

#ifdef __cplusplus
}
#endif

#endif /* AMR_EKF_H */
