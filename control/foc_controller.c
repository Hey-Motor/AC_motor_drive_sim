#include "foc_controller.h"
#include "transforms.h"
#include <math.h>
extern double g_debug_time;

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
    //double id_ref = 0.0;
    /* 为 PEBO-DREM 提供额外电气激励：
    * 在 d 轴注入一个慢变小幅正弦
    */
    double id_ref = 0;//1.0 * sin(2.0 * 3.14159265358979323846 * 3.0 * g_debug_time);
    double iq_ref;
    double Ts = p->sim.Ts_ctrl;

    park(adc->i_alpha, adc->i_beta, theta_e_fb, &id, &iq);

    {
        double err_spd = speed_ref - omega_m_fb;
        s->spd_int += p->ctrl.ki_spd * Ts * err_spd;
        s->spd_int = clamp(s->spd_int, -p->ctrl.iq_limit, p->ctrl.iq_limit);

        //nlrs那个仿真必须注入q轴电流正弦扰动  Rs才能收敛 
        //若w_ref=30rad/s，则注入幅值2.0；  100rad/s 注入0.5

        //2011-inoue 论文可以不注入（但是必须加载）如果空载必须注入才能电阻收敛
        iq_ref = p->ctrl.kp_spd * err_spd + s->spd_int;// + 2.0 * sin(2.0 * 3.14159265358979323846 * 45.0 * g_debug_time);
        iq_ref = clamp(iq_ref, -p->ctrl.iq_limit, p->ctrl.iq_limit);
    }

    {
        double err_d = id_ref - id;
        double err_q = iq_ref - iq;
        double omega_e = p->motor.pole_pairs * omega_m_fb;
        double vd, vq;

        s->id_int += p->ctrl.ki_id * Ts * err_d;
        s->iq_int += p->ctrl.ki_iq * Ts * err_q;

        s->id_int = clamp(s->id_int, -p->ctrl.vdq_limit, p->ctrl.vdq_limit);
        s->iq_int = clamp(s->iq_int, -p->ctrl.vdq_limit, p->ctrl.vdq_limit);

        vd = p->ctrl.kp_id * err_d + s->id_int - omega_e * p->motor.Lq * iq;
        vq = p->ctrl.kp_iq * err_q + s->iq_int + omega_e * (p->motor.Ld * id + p->motor.psi_f);

        circle_limit(&vd, &vq, p->ctrl.vdq_limit);
        inv_park(vd, vq, theta_e_fb, u_alpha_cmd, u_beta_cmd);
    }

    if (id_ref_out) *id_ref_out = id_ref;
    if (iq_ref_out) *iq_ref_out = iq_ref;
}