#include "foc_controller.h"
#include "../config/config.h"
#include "transforms.h"
#include <math.h>
extern double g_debug_time;
static int g_speed_gain_sched_enabled = 0;

static double clamp(double x, double xmin, double xmax)
{
    if (x > xmax) return xmax;
    if (x < xmin) return xmin;
    return x;
}

void foc_init(FocState *s)
{
    s->id_int  = 0.0;
    s->iq_int  = 0.0;
    s->spd_int = 0.0;
    s->flux_int = 0.0;
}

void foc_set_speed_gain_sched_enabled(int enabled)
{
    g_speed_gain_sched_enabled = enabled ? 1 : 0;
}

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
)
{
    double id, iq;
    double id_ref = 0.0;
    double iq_ref;
    double Ts = p->sim.Ts_ctrl;

    park(adc->i_alpha, adc->i_beta, theta_e_fb, &id, &iq);

    {
        double err_spd = speed_ref - omega_m_fb;
#if (MOTOR_TYPE_DEFAULT == MOTOR_TYPE_IM)
        double omega_e = p->motor_im.pole_pairs * omega_m_fb;
        double Lsigma = p->motor_im.Ls - (p->motor_im.Lm * p->motor_im.Lm) / p->motor_im.Lr;
        double err_flux = p->ctrl_im.psi_r_ref - id * p->motor_im.Lm;
        double kp_spd_use = p->ctrl_im.kp_spd;
        double ki_spd_use = p->ctrl_im.ki_spd;
        double vd, vq;
        double err_d, err_q;

        if (g_speed_gain_sched_enabled && fabs(omega_m_fb) > 90.0) {
            kp_spd_use = 0.6 * p->ctrl_im.kp_spd;
            ki_spd_use = 0.6 * p->ctrl_im.ki_spd;
        }

        s->spd_int += ki_spd_use * Ts * err_spd;
        s->spd_int = clamp(s->spd_int, -p->ctrl_im.iq_limit, p->ctrl_im.iq_limit);
        iq_ref = kp_spd_use * err_spd + s->spd_int;
        iq_ref = clamp(iq_ref, -p->ctrl_im.iq_limit, p->ctrl_im.iq_limit);

        s->flux_int += p->ctrl_im.ki_flux * Ts * err_flux;
        s->flux_int = clamp(s->flux_int, -p->ctrl_im.id_limit, p->ctrl_im.id_limit);
        id_ref = (p->ctrl_im.kp_flux * err_flux + s->flux_int) / (p->motor_im.Lm + 1e-6);
        id_ref = clamp(id_ref, -p->ctrl_im.id_limit, p->ctrl_im.id_limit);

        err_d = id_ref - id;
        err_q = iq_ref - iq;

        s->id_int += p->ctrl_im.ki_id * Ts * err_d;
        s->iq_int += p->ctrl_im.ki_iq * Ts * err_q;
        s->id_int = clamp(s->id_int, -p->ctrl_im.vdq_limit, p->ctrl_im.vdq_limit);
        s->iq_int = clamp(s->iq_int, -p->ctrl_im.vdq_limit, p->ctrl_im.vdq_limit);

        vd = p->ctrl_im.kp_id * err_d + s->id_int - omega_e * Lsigma * iq;
        vq = p->ctrl_im.kp_iq * err_q + s->iq_int + omega_e * Lsigma * id;

        circle_limit(&vd, &vq, p->ctrl_im.vdq_limit);
        inv_park(vd, vq, theta_e_fb, u_alpha_cmd, u_beta_cmd);
#else
        {
            double err_d, err_q;
            double omega_e = p->motor.pole_pairs * omega_m_fb;
            double vd, vq;

            s->spd_int += p->ctrl.ki_spd * Ts * err_spd;
            s->spd_int = clamp(s->spd_int, -p->ctrl.iq_limit, p->ctrl.iq_limit);
            iq_ref = p->ctrl.kp_spd * err_spd + s->spd_int;
            iq_ref = clamp(iq_ref, -p->ctrl.iq_limit, p->ctrl.iq_limit);

            err_d = id_ref - id;
            err_q = iq_ref - iq;

            s->id_int += p->ctrl.ki_id * Ts * err_d;
            s->iq_int += p->ctrl.ki_iq * Ts * err_q;

            s->id_int = clamp(s->id_int, -p->ctrl.vdq_limit, p->ctrl.vdq_limit);
            s->iq_int = clamp(s->iq_int, -p->ctrl.vdq_limit, p->ctrl.vdq_limit);

            vd = p->ctrl.kp_id * err_d + s->id_int - omega_e * p->motor.Lq * iq;
            vq = p->ctrl.kp_iq * err_q + s->iq_int + omega_e * (p->motor.Ld * id + p->motor.psi_f);

            circle_limit(&vd, &vq, p->ctrl.vdq_limit);
            inv_park(vd, vq, theta_e_fb, u_alpha_cmd, u_beta_cmd);
        }
#endif
    }

    if (id_ref_out) *id_ref_out = id_ref;
    if (iq_ref_out) *iq_ref_out = iq_ref;
}
