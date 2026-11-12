/*
 * Cinemática de bases holonómicas de 4 ruedas. Misma convención que
 * AMR_LLC/ros2/amr_mqtt_bridge/amr_mqtt_bridge/omni.py.
 *
 * Mecanum (el robot real): ruedas rectas en (+-lx, +-ly), rodillos a 45 deg en X
 * vistos desde arriba, w_i positiva = la rueda rueda hacia +x:
 *
 *   fl:  w r = vx - vy - (lx + ly) wz        rr:  w r = vx - vy + (lx + ly) wz
 *   rl:  w r = vx + vy - (lx + ly) wz        fr:  w r = vx + vy + (lx + ly) wz
 *
 * Omni en X (amr_omni_init): rueda i en (L cos a_i, L sin a_i), eje radial,
 * tracción tangencial: v_i = -sin(a_i) vx + cos(a_i) vy + L wz,  w_i = v_i / r.
 *
 * Orden de ruedas del contrato: fl, rl, rr, fr (ruedas 1, 2, 3, 4).
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
	float L;                       /* brazo de palanca de wz [m] */
	float J[AMR_OMNI_WHEELS][3];   /* v_superficie = J [vx vy wz] */
	float Jp[3][AMR_OMNI_WHEELS];  /* pseudo-inversa (J^T J)^-1 J^T */
} amr_omni_t;

/* angles_rad: ángulo de cada rueda alrededor del centro. false si la
 * configuración es singular o los parámetros no son positivos. */
bool amr_omni_init(amr_omni_t *k, float wheel_radius, float center_distance,
		   const float angles_rad[AMR_OMNI_WHEELS]);

/* Mecanum con ruedas en (+-half_length, +-half_width). false si algún
 * parámetro no es positivo. */
bool amr_mecanum_init(amr_omni_t *k, float wheel_radius, float half_length, float half_width);

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
