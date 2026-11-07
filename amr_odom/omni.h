/*
 * Cinemática de la base omnidireccional de 4 ruedas (configuración X a 45 deg).
 * Misma convención que AMR_LLC/ros2/amr_mqtt_bridge/amr_mqtt_bridge/omni.py:
 *
 *   rueda i en (L cos a_i, L sin a_i), eje radial, tracción tangencial
 *   u_i = (-sin a_i, cos a_i):    v_i = -sin(a_i) vx + cos(a_i) vy + L wz,  w_i = v_i / r
 *
 * Orden de ruedas del contrato: fl (45), rl (135), rr (225), fr (315).
 * C99 portable, sin RTOS.
 */
#ifndef AMR_OMNI_H
#define AMR_OMNI_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AMR_OMNI_WHEELS 4

typedef struct {
	float r;                       /* radio de rueda [m] */
	float L;                       /* centro -> rueda [m] */
	float J[AMR_OMNI_WHEELS][3];   /* v_superficie = J [vx vy wz] */
	float Jp[3][AMR_OMNI_WHEELS];  /* pseudo-inversa (J^T J)^-1 J^T */
} amr_omni_t;

/* angles_rad: ángulo de cada rueda alrededor del centro. false si la
 * configuración es singular o los parámetros no son positivos. */
bool amr_omni_init(amr_omni_t *k, float wheel_radius, float center_distance,
		   const float angles_rad[AMR_OMNI_WHEELS]);

/* Velocidad del robot (base_footprint) -> velocidad angular de cada rueda [rad/s]. */
void amr_omni_inverse(const amr_omni_t *k, float vx, float vy, float wz,
		      float wheel_rad_s[AMR_OMNI_WHEELS]);

/* Velocidades (o desplazamientos) angulares de rueda -> (vx, vy, wz) por mínimos cuadrados. */
void amr_omni_forward(const amr_omni_t *k, const float wheel_rad_s[AMR_OMNI_WHEELS],
		      float *vx, float *vy, float *wz);

/* Error RMS de la directa: mide qué tan inconsistentes son las 4 ruedas
 * (deslizamiento). 0 para mediciones cinemáticamente perfectas. [m/s] */
float amr_omni_residual(const amr_omni_t *k, const float wheel_rad_s[AMR_OMNI_WHEELS]);

#ifdef __cplusplus
}
#endif

#endif /* AMR_OMNI_H */
