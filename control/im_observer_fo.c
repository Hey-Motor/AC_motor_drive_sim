#include "im_observer_fo.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double g_i_alpha_hat = 0.0;
static double g_i_beta_hat = 0.0;
static double g_psi_r_alpha_hat = 0.8;
static double g_psi_r_beta_hat = 0.0;
static double g_wr_hat = 0.0;
static double g_wr_int = 0.0;
static double g_omega_adapt_err_f = 0.0;
static double g_rs_hat = 1.5;
static double g_rs_int = 0.0;
static double g_kp_rs = 0.05;
static double g_ki_rs = 5.0;
static double g_rs_err_f = 0.0;
static double g_rs_nominal = 0.567;
static int g_use_real_speed = 0;
static int g_adapt_gain_sched_enabled = 0;
static int g_adapt_err_lpf_enabled = 0;
static int g_rs_adapt_enabled = 1;

typedef struct {
    double di_a_hat;
    double di_b_hat;
    double dpsi_ra_hat;
    double dpsi_rb_hat;
} ObsDeriv;

static double wrap_pm_pi(double x)
{
    while (x > M_PI) x -= 2.0 * M_PI;
    while (x < -M_PI) x += 2.0 * M_PI;
    return x;
}

void im_observer_fo_init(ObserverState *s)
{
    if (s != NULL) {
        memset(s, 0, sizeof(*s));
        s->initialized = 1;
    }

    g_i_alpha_hat = 0.0;
    g_i_beta_hat = 0.0;
    g_psi_r_alpha_hat = 0.8;
    g_psi_r_beta_hat = 0.0;
    g_wr_hat = 0.0;
    g_wr_int = 0.0;
    g_omega_adapt_err_f = 0.0;
    g_rs_hat = 1.5;
    g_rs_int = 0.0;
    g_rs_err_f = 0.0;
    g_rs_nominal = 0.567;
}

void im_observer_fo_set_speed_init(double omega_e_init)
{
    g_wr_hat = omega_e_init;
    g_wr_int = omega_e_init;
}

void im_observer_fo_set_use_real_speed(int use_real_speed)
{
    g_use_real_speed = use_real_speed ? 1 : 0;
}

void im_observer_fo_set_adapt_gain_sched_enabled(int enabled)
{
    g_adapt_gain_sched_enabled = enabled ? 1 : 0;
}

void im_observer_fo_set_adapt_err_lpf_enabled(int enabled)
{
    g_adapt_err_lpf_enabled = enabled ? 1 : 0;
}

void im_observer_fo_set_rs_adapt_enabled(int enabled)
{
    g_rs_adapt_enabled = enabled ? 1 : 0;
}

void im_observer_fo_set_rs_adapt_gains(double kp_r, double ki_r)
{
    g_kp_rs = kp_r;
    g_ki_rs = ki_r;
}

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
)
{
    double Tr;
    double sigma;
    double wr_model;
    double rs_used;
    double a, b, c, d;
    double ei_a, ei_b;
    double J_psira_hat, J_psirb_hat;
    double gi = 300.0;
    double omega_adapt_err;
    double psi_mag;
    double kp_w, ki_w;
    double rs_adapt_err;

    if ((s == NULL) || (out == NULL) || (mp == NULL)) {
        return;
    }
    if (!s->initialized) {
        im_observer_fo_init(s);
    }

    Tr = mp->Lr / mp->Rr;
    sigma = 1.0 - (mp->Lm * mp->Lm) / (mp->Ls * mp->Lr);
    wr_model = g_use_real_speed ? (mp->pole_pairs * omega_m_real) : g_wr_hat;

    rs_used = g_rs_hat;
    if (rs_used < 0.05) {
        rs_used = 0.05;
    } else if (rs_used > 3.0) {
        rs_used = 3.0;
    }
    a = -(rs_used / (sigma * mp->Ls) + (mp->Lm * mp->Lm) / (sigma * mp->Ls * mp->Lr * Tr));
    b = 1.0 / (sigma * mp->Ls);
    c = mp->Lm / (sigma * mp->Ls * mp->Lr * Tr);
    d = mp->Lm / Tr;

    ei_a = i_alpha - g_i_alpha_hat;
    ei_b = i_beta - g_i_beta_hat;

    J_psira_hat = -g_psi_r_beta_hat;
    J_psirb_hat = g_psi_r_alpha_hat;

    {
        ObsDeriv k1, k2, k3, k4;
        double ia2, ib2, psia2, psib2;
        double ia3, ib3, psia3, psib3;
        double ia4, ib4, psia4, psib4;
        double ei_a2, ei_b2, ei_a3, ei_b3, ei_a4, ei_b4;
        double J2_a, J2_b, J3_a, J3_b, J4_a, J4_b;

        k1.di_a_hat = a * g_i_alpha_hat + c * g_psi_r_alpha_hat - c * Tr * wr_model * J_psira_hat + b * u_alpha
                    + gi * ei_a;
        k1.di_b_hat = a * g_i_beta_hat  + c * g_psi_r_beta_hat  - c * Tr * wr_model * J_psirb_hat + b * u_beta
                    + gi * ei_b;
        k1.dpsi_ra_hat = d * i_alpha - (1.0 / Tr) * g_psi_r_alpha_hat + wr_model * J_psira_hat;
        k1.dpsi_rb_hat = d * i_beta  - (1.0 / Tr) * g_psi_r_beta_hat  + wr_model * J_psirb_hat;

        ia2 = g_i_alpha_hat + 0.5 * Ts * k1.di_a_hat;
        ib2 = g_i_beta_hat + 0.5 * Ts * k1.di_b_hat;
        psia2 = g_psi_r_alpha_hat + 0.5 * Ts * k1.dpsi_ra_hat;
        psib2 = g_psi_r_beta_hat + 0.5 * Ts * k1.dpsi_rb_hat;
        ei_a2 = i_alpha - ia2;
        ei_b2 = i_beta - ib2;
        J2_a = -psib2;
        J2_b = psia2;
        k2.di_a_hat = a * ia2 + c * psia2 - c * Tr * wr_model * J2_a + b * u_alpha + gi * ei_a2;
        k2.di_b_hat = a * ib2 + c * psib2 - c * Tr * wr_model * J2_b + b * u_beta + gi * ei_b2;
        k2.dpsi_ra_hat = d * i_alpha - (1.0 / Tr) * psia2 + wr_model * J2_a;
        k2.dpsi_rb_hat = d * i_beta  - (1.0 / Tr) * psib2 + wr_model * J2_b;

        ia3 = g_i_alpha_hat + 0.5 * Ts * k2.di_a_hat;
        ib3 = g_i_beta_hat + 0.5 * Ts * k2.di_b_hat;
        psia3 = g_psi_r_alpha_hat + 0.5 * Ts * k2.dpsi_ra_hat;
        psib3 = g_psi_r_beta_hat + 0.5 * Ts * k2.dpsi_rb_hat;
        ei_a3 = i_alpha - ia3;
        ei_b3 = i_beta - ib3;
        J3_a = -psib3;
        J3_b = psia3;
        k3.di_a_hat = a * ia3 + c * psia3 - c * Tr * wr_model * J3_a + b * u_alpha + gi * ei_a3;
        k3.di_b_hat = a * ib3 + c * psib3 - c * Tr * wr_model * J3_b + b * u_beta + gi * ei_b3;
        k3.dpsi_ra_hat = d * i_alpha - (1.0 / Tr) * psia3 + wr_model * J3_a;
        k3.dpsi_rb_hat = d * i_beta  - (1.0 / Tr) * psib3 + wr_model * J3_b;

        ia4 = g_i_alpha_hat + Ts * k3.di_a_hat;
        ib4 = g_i_beta_hat + Ts * k3.di_b_hat;
        psia4 = g_psi_r_alpha_hat + Ts * k3.dpsi_ra_hat;
        psib4 = g_psi_r_beta_hat + Ts * k3.dpsi_rb_hat;
        ei_a4 = i_alpha - ia4;
        ei_b4 = i_beta - ib4;
        J4_a = -psib4;
        J4_b = psia4;
        k4.di_a_hat = a * ia4 + c * psia4 - c * Tr * wr_model * J4_a + b * u_alpha + gi * ei_a4;
        k4.di_b_hat = a * ib4 + c * psib4 - c * Tr * wr_model * J4_b + b * u_beta + gi * ei_b4;
        k4.dpsi_ra_hat = d * i_alpha - (1.0 / Tr) * psia4 + wr_model * J4_a;
        k4.dpsi_rb_hat = d * i_beta  - (1.0 / Tr) * psib4 + wr_model * J4_b;

        g_i_alpha_hat += Ts * (k1.di_a_hat + 2.0 * k2.di_a_hat + 2.0 * k3.di_a_hat + k4.di_a_hat) / 6.0;
        g_i_beta_hat += Ts * (k1.di_b_hat + 2.0 * k2.di_b_hat + 2.0 * k3.di_b_hat + k4.di_b_hat) / 6.0;
        g_psi_r_alpha_hat += Ts * (k1.dpsi_ra_hat + 2.0 * k2.dpsi_ra_hat + 2.0 * k3.dpsi_ra_hat + k4.dpsi_ra_hat) / 6.0;
        g_psi_r_beta_hat += Ts * (k1.dpsi_rb_hat + 2.0 * k2.dpsi_rb_hat + 2.0 * k3.dpsi_rb_hat + k4.dpsi_rb_hat) / 6.0;
    }

    /* ========================= 速度自适应律 =========================
     * 误差信号：e_w = e_i^T * J * psi_hat / |psi_hat|^2
     * 估计律：   w_hat = -(Kp + Ki/s) * e_w
     * 说明：
     * 1) e_i = i - i_hat
     * 2) J*psi_hat 相当于把磁链向量旋转 90°，用于构造速度相关误差
     * =============================================================== */
    psi_mag = hypot(g_psi_r_alpha_hat, g_psi_r_beta_hat);
    omega_adapt_err = (ei_a * (-g_psi_r_beta_hat) + ei_b * g_psi_r_alpha_hat) / (psi_mag * psi_mag + 1e-4);
    if (g_adapt_err_lpf_enabled) {
        const double tau = 0.02;
        const double alpha = Ts / (tau + Ts);
        g_omega_adapt_err_f += alpha * (omega_adapt_err - g_omega_adapt_err_f);
        omega_adapt_err = g_omega_adapt_err_f;
    }

    if (g_adapt_gain_sched_enabled) {
        if (fabs(wr_model / mp->pole_pairs) <= 90.0) {
            kp_w = 20.0;
            ki_w = 20.0;
        } else {
            kp_w = 6.0;
            ki_w = 8.0;
        }
    } else {
        kp_w = 20.0;
        ki_w = 20.0;
    }

    g_wr_int += Ts * (-ki_w * omega_adapt_err);
    {
        double wr_pi = -kp_w * omega_adapt_err + g_wr_int;
        g_wr_hat = wr_pi;
    }

    {   /* ========================= Rs自适应误差 =========================
         * e_R = (e_i^T * i_hat) / |i_hat|^2
         * 先做归一化，再一阶低通，降低数值抖动。
         * =============================================================== */
        double reg_den = g_i_alpha_hat * g_i_alpha_hat + g_i_beta_hat * g_i_beta_hat + 1e-4;
        const double tau_r = 0.02;
        const double alpha_r = Ts / (tau_r + Ts);
        rs_adapt_err = (ei_a * g_i_alpha_hat + ei_b * g_i_beta_hat) / reg_den;
        g_rs_err_f += alpha_r * (rs_adapt_err - g_rs_err_f);
    }
    if (g_rs_adapt_enabled) {
        /* ========================= Rs自适应律 =========================
         * R_hat = -(Kp + Ki/s) * e_R
         * 这里加了一个很小的泄漏项到名义值，防止纯积分长期漂移：
         *     d(int_R)/dt = -Ki*e_R - k_leak*(R_hat - R_nominal)
         * 最后再做投影限幅，保证估计电阻保持在物理范围内。
         * =============================================================== */
        const double k_leak = 1.0;
        g_rs_int += Ts * (-g_ki_rs * g_rs_err_f - k_leak * (g_rs_hat - g_rs_nominal));
        g_rs_hat = -g_kp_rs * g_rs_err_f + g_rs_int;
        if (g_rs_hat < 0.05) {
            g_rs_hat = 0.05;
            if (g_rs_int < 0.05) g_rs_int = 0.05;
        } else if (g_rs_hat > 3.0) {
            g_rs_hat = 3.0;
            if (g_rs_int > 3.0) g_rs_int = 3.0;
        }
    }

    memset(out, 0, sizeof(*out));
    out->psi_alpha_hat = g_psi_r_alpha_hat;
    out->psi_beta_hat = g_psi_r_beta_hat;
    out->phi_hat = hypot(g_psi_r_alpha_hat, g_psi_r_beta_hat);

    out->theta_e_hat = wrap_pm_pi(atan2(g_psi_r_beta_hat, g_psi_r_alpha_hat));
    out->omega_e_hat = g_wr_hat;

    out->z21 = ei_a;
    out->z22 = ei_b;
    out->yreg = omega_adapt_err;
    out->R_hat = g_rs_hat;
    out->eta1_hat = g_rs_err_f;
    out->eta2_hat = g_kp_rs;
    out->beta_hat = g_ki_rs;

    s->psi_alpha_hat = out->psi_alpha_hat;
    s->psi_beta_hat = out->psi_beta_hat;
    s->phi_hat = out->phi_hat;
    s->theta_prev = out->theta_e_hat;
    s->omega_e_hat_f = out->omega_e_hat;
}
