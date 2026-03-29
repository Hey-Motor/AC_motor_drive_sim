#include "observer_piippo2008.h"

#include <math.h>
#include <string.h>

#include "../config/config.h"

/*
 * Piippo / Hinkkanen / Luomi 2008 adaptive observer (baseline)
 * -------------------------------------------------------------
 * 本阶段只实现论文中的基础 speed-adaptive observer：式 (4)~(10)。
 * 不实现 HF signal injection（式 (17)~(22)），后续再补。
 *
 * 另外按本轮要求：主控制 FOC 继续由传感器/真值驱动，
 * 本 observer 仅用于独立评估收敛表现。
 */

#if (USE_IPMSM == 1)
#define PIIPPO_R_S       (1.025)
#define PIIPPO_L_D       (9.0e-3)
#define PIIPPO_L_Q       (18.5e-3)
#define PIIPPO_PSI_PM    (0.249)
#else
#define PIIPPO_R_S       (1.025)
#define PIIPPO_L_D       (9.0e-3)
#define PIIPPO_L_Q       (9.0e-3)
#define PIIPPO_PSI_PM    (0.249)
#endif

#define PIIPPO_POLE_PAIRS        (3.0)
#define PIIPPO_OMEGA_BASE_E      (471.23889803846896) /* 2*pi*75 rad/s */

/* 速度自适应带宽 alpha_o (rad/s)，对应式(10) */
#define PIIPPO_ALPHA_O           (120.0)
#define PIIPPO_KP                (2.0 * PIIPPO_ALPHA_O / PIIPPO_PSI_PM)
#define PIIPPO_KI                ((PIIPPO_ALPHA_O * PIIPPO_ALPHA_O) / PIIPPO_PSI_PM)

#define PIIPPO_ADAPT_SIGN        (-1.0)  /* 可调：+1/-1 用于校正自适应律符号 */
#define PIIPPO_FE_SIGN           (1.0)   /* 可调：+1/-1 用于 F_eps 方向 */

/* 观测器增益：支持常值增益和 speed-dependent gain */
#define PIIPPO_GAIN_CONST        (0)
#define PIIPPO_GAIN_SPEED_DEP    (1)
#define PIIPPO_GAIN_MODE         PIIPPO_GAIN_SPEED_DEP

/* 常值增益（论文讨论值之一） */
#define PIIPPO_LAMBDA1_CONST     (-0.5 * PIIPPO_R_S)
#define PIIPPO_LAMBDA2_CONST     (0.0)

/* speed-dependent 增益，式 (16a)(16b)
 * 论文中的 λ' 是独立设计常数，这里显式设为常数（单位等价于欧姆）。
 * 默认值 2.05 仅是工程初值，数值上约等于 2*Rs(1.025)，但不再与 Rs 绑定。 */
#define PIIPPO_LAMBDA_PRIME      (2.05)
#define PIIPPO_OMEGA_LAMBDA      (1.0 * PIIPPO_OMEGA_BASE_E)

/* 数值保护 */
#define PIIPPO_PI                (3.14159265358979323846)
#define PIIPPO_TWO_PI            (6.28318530717958647692)
#define PIIPPO_OMEGA_E_MAX       (2000.0)
#define PIIPPO_INT_MAX           (2000.0)

static double g_theta_e_hat = 0.0;
static double g_omega_e_hat = 0.0;
static double g_speed_int = 0.0;

/* 式(4) 中状态变量：估计定子磁链（估计转子坐标系） */
static double g_psi_d_hat = PIIPPO_PSI_PM;
static double g_psi_q_hat = 0.0;

static double wrap_pm_pi(double x)
{
    while (x > PIIPPO_PI) {
        x -= PIIPPO_TWO_PI;
    }
    while (x < -PIIPPO_PI) {
        x += PIIPPO_TWO_PI;
    }
    return x;
}

static double clampd(double x, double xmin, double xmax)
{
    if (x < xmin) {
        return xmin;
    }
    if (x > xmax) {
        return xmax;
    }
    return x;
}

static void select_observer_gain(double omega_e_hat, double *lambda1, double *lambda2)
{
#if (PIIPPO_GAIN_MODE == PIIPPO_GAIN_CONST)
    *lambda1 = PIIPPO_LAMBDA1_CONST;
    *lambda2 = PIIPPO_LAMBDA2_CONST;
#else
    double abs_pu = fabs(omega_e_hat) / PIIPPO_OMEGA_LAMBDA;

    if (abs_pu <= 1.0) {
        *lambda1 = PIIPPO_LAMBDA_PRIME * abs_pu;
        *lambda2 = PIIPPO_LAMBDA_PRIME * (omega_e_hat / PIIPPO_OMEGA_LAMBDA);
    } else {
        *lambda1 = PIIPPO_LAMBDA_PRIME;
        *lambda2 = PIIPPO_LAMBDA_PRIME * ((omega_e_hat >= 0.0) ? 1.0 : -1.0);
    }
#endif
}

void observer_piippo2008_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}

void observer_piippo2008_init(ObserverState *s)
{
    if (s != NULL) {
        memset(s, 0, sizeof(*s));
        s->initialized = 1;
    }

    g_theta_e_hat = 0.0;
    g_omega_e_hat = 0.0;
    g_speed_int = 0.0;

    g_psi_d_hat = PIIPPO_PSI_PM;
    g_psi_q_hat = 0.0;
}

void observer_piippo2008_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    double c, sn;
    double i_d_ref, i_q_ref;
    double u_d_ref, u_q_ref;

    double i_d_hat, i_q_hat;
    double i_tilde_d, i_tilde_q;

    double lambda1, lambda2;
    double F_eps;
    double dpsi_d, dpsi_q;

    if ((s == NULL) || (out == NULL)) {
        return;
    }
    if (!s->initialized) {
        observer_piippo2008_init(s);
    }

    c = cos(g_theta_e_hat);
    sn = sin(g_theta_e_hat);

    /* measured quantities in estimated rotor frame (论文 i'_s, u'_s) */
    i_d_ref = c * i_alpha + sn * i_beta;
    i_q_ref = -sn * i_alpha + c * i_beta;

    u_d_ref = c * u_alpha + sn * u_beta;
    u_q_ref = -sn * u_alpha + c * u_beta;

    /* 式(5): current estimate from flux estimate */
    i_d_hat = (g_psi_d_hat - PIIPPO_PSI_PM) / PIIPPO_L_D;
    i_q_hat = g_psi_q_hat / PIIPPO_L_Q;

    /* 式(6): current estimation error */
    i_tilde_d = i_d_ref - i_d_hat;
    i_tilde_q = i_q_ref - i_q_hat;

    /* 式(7) 与式(16): lambda = lambda1*I + lambda2*J */
    select_observer_gain(g_omega_e_hat, &lambda1, &lambda2);

    /* 式(8): adaptation error, C1=[0 Lq] */
    F_eps = PIIPPO_FE_SIGN * PIIPPO_L_Q * i_tilde_q;

    /* 式(9): PI speed adaptation */
    g_speed_int += Ts * (PIIPPO_ADAPT_SIGN * PIIPPO_KI * F_eps);
    g_speed_int = clampd(g_speed_int, -PIIPPO_INT_MAX, PIIPPO_INT_MAX);

    g_omega_e_hat = PIIPPO_ADAPT_SIGN * PIIPPO_KP * F_eps + g_speed_int;
    g_omega_e_hat = clampd(g_omega_e_hat, -PIIPPO_OMEGA_E_MAX, PIIPPO_OMEGA_E_MAX);

    g_theta_e_hat = wrap_pm_pi(g_theta_e_hat + Ts * g_omega_e_hat);

    /* 式(4): observer model (estimated rotor frame) */
    dpsi_d = u_d_ref
           - PIIPPO_R_S * i_d_hat
           + g_omega_e_hat * g_psi_q_hat
           + lambda1 * i_tilde_d
           - lambda2 * i_tilde_q;

    dpsi_q = u_q_ref
           - PIIPPO_R_S * i_q_hat
           - g_omega_e_hat * g_psi_d_hat
           + lambda2 * i_tilde_d
           + lambda1 * i_tilde_q;

    g_psi_d_hat += Ts * dpsi_d;
    g_psi_q_hat += Ts * dpsi_q;

    memset(out, 0, sizeof(*out));
    out->theta_e_hat = g_theta_e_hat;
    out->omega_e_hat = g_omega_e_hat;

    /* 仅用于调试输出 */
    out->psi_alpha_hat = c * g_psi_d_hat - sn * g_psi_q_hat;
    out->psi_beta_hat = sn * g_psi_d_hat + c * g_psi_q_hat;
    out->phi_hat = hypot(out->psi_alpha_hat, out->psi_beta_hat);

    out->z21 = i_tilde_d;
    out->z22 = i_tilde_q;
    out->xi1 = lambda1;
    out->xi2 = lambda2;
    out->yreg = F_eps;

    s->psi_alpha_hat = out->psi_alpha_hat;
    s->psi_beta_hat = out->psi_beta_hat;
    s->phi_hat = out->phi_hat;
    s->theta_prev = g_theta_e_hat;
    s->omega_e_hat_f = g_omega_e_hat;
    s->pll_theta = g_theta_e_hat;
    s->pll_omega_i = g_speed_int;
}
