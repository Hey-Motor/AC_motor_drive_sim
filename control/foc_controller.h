#ifndef FOC_CONTROLLER_H
#define FOC_CONTROLLER_H

#include "../core/types.h"

void foc_init(FocState *s);

void foc_step(
    double speed_ref,
    double theta_e_fb,
    double omega_m_fb,
    const AdcSample *adc,
    const Params *p,
    FocState *s,
    double *u_alpha_cmd,
    double *u_beta_cmd,
    double *id_ref_out,
    double *iq_ref_out
);
extern double g_debug_time;
#endif