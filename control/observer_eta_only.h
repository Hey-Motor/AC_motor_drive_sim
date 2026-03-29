#ifndef OBSERVER_ETA_ONLY_H
#define OBSERVER_ETA_ONLY_H

#include "../core/types.h"

void observer_eta_only_init(ObserverState *s);
void observer_eta_only_set_debug_truth(double theta_e_true, double omega_e_true);
void observer_eta_only_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif
