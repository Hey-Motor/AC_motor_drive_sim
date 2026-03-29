#ifndef OBSERVER_IPMSM_NLRS_H
#define OBSERVER_IPMSM_NLRS_H

#include "observer_select.h"

void observer_ipmsm_nlrs_init(ObserverState *s);
void observer_ipmsm_nlrs_set_debug_truth(double theta_e_true, double omega_e_true);

void observer_ipmsm_nlrs_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif
