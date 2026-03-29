#ifndef OBSERVER_PEBO_DREM_H
#define OBSERVER_PEBO_DREM_H

#include "observer_select.h"

/* Bazylev/Pyrkin/Bobtsov 2018:
 * PMSM with unknown stator resistance
 * PEBO + DREM + hybrid theta0 reconstruction
 */

void observer_pebo_drem_init(ObserverState *s);
void observer_pebo_drem_set_debug_truth(double theta_e_true, double omega_e_true);

void observer_pebo_drem_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif