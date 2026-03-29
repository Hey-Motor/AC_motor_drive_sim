#ifndef OBSERVER_PIIPPO2008_H
#define OBSERVER_PIIPPO2008_H

#include "../core/types.h"

void observer_piippo2008_init(ObserverState *s);
void observer_piippo2008_set_debug_truth(double theta_e_true, double omega_e_true);
void observer_piippo2008_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif
