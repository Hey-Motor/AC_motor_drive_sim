#ifndef IM_OBSERVER_FO_H
#define IM_OBSERVER_FO_H

#include "../core/types.h"

void im_observer_fo_init(ObserverState *s);
void im_observer_fo_set_speed_init(double omega_e_init);
void im_observer_fo_set_use_real_speed(int use_real_speed);
void im_observer_fo_set_adapt_gain_sched_enabled(int enabled);
void im_observer_fo_set_adapt_err_lpf_enabled(int enabled);
void im_observer_fo_set_rs_adapt_enabled(int enabled);
void im_observer_fo_set_rs_adapt_gains(double kp_r, double ki_r);
void im_observer_fo_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    const ImMotorParams *mp,
    double omega_m_real,
    ObserverOutput *out
);

#endif
