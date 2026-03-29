#ifndef OBSERVER_BOBTSOV2015_H
#define OBSERVER_BOBTSOV2015_H

#include "observer_select.h"

void observer_bobtsov2015_init(ObserverState *s);
void observer_bobtsov2015_set_debug_truth(double theta_e_true, double omega_e_true);

void observer_bobtsov2015_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif