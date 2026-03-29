#include "observer_select.h"
#include "../config/config.h"
#include <math.h>

/* =========================================================
 * 非线性磁链观测器 + 在线永磁磁链估计
 *
 * 核心连续时间方程：
 *   dot(Psi_hat) = u - R i - 2 q (Psi_hat - L i) (|Psi_hat - L i|^2 - Phi_hat^2)
 *   dot(Phi_hat) = q Phi_hat (|Psi_hat - L i|^2 - Phi_hat^2)
 *
 * 其中：
 *   Psi_hat = 总磁链估计
 *   Phi_hat = 永磁磁链幅值估计
 *   z = Psi_hat - L i
 * z 的方向就是转子磁链方向，因此可以从 z 提取电角度。
 *
 * 当前文件做了两种“速度提取”方式：
 *
 * 1) OBS_SPEED_ATAN_DIFF
 *    theta_hat = atan2(z_beta, z_alpha)
 *    omega_hat = d(theta_hat)/dt，再加一阶低通
 *
 * 2) OBS_SPEED_PLL
 *    用磁链向量 z 进入正交 PLL
 *    PLL 直接输出 theta_hat, omega_hat
 *
 * 推荐：
 *   先用 ATAN_DIFF 验证趋势，再切 PLL 做正式结果
 * ========================================================= */


/* =========================
 * 观测器参数
 * 当前先写在本文件，方便快速调参
 * 后面若你需要，可以再移到 params.c
 * ========================= */
#define OBS_R               (1.025)
#define OBS_L               (9e-3)      /* SPMSM 用标量 L */
#define OBS_Q               (2000.0)

#define OBS_PHI_INIT        (0.20)
#define OBS_PHI_MIN         (1e-4)

/* atan2+微分法的一阶低通带宽 */
#define OBS_OMEGA_LPF_WC    (400.0)

/* PLL 参数
 * 建议先从这个量级开始
 * 不合适再调
 */
#define OBS_PLL_KP          (300.0)
#define OBS_PLL_KI          (30000.0)

/* 误差归一化防止幅值过小时数值发散 */
#define OBS_Z_NORM_EPS      (1e-6)


static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < -3.14159265358979323846) x += 2.0 * 3.14159265358979323846;
    return x;
}


/* 旧接口保留，当前真实观测器里不用它 */
void observer_fluxphi_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}


void observer_fluxphi_init(ObserverState *s)
{
    s->psi_alpha_hat = OBS_PHI_INIT;
    s->psi_beta_hat  = 0.0;
    s->phi_hat       = OBS_PHI_INIT;

    s->theta_prev    = 0.0;
    s->omega_e_hat_f = 0.0;

    s->pll_theta     = 0.0;
    s->pll_omega_i   = 0.0;

    s->initialized   = 1;
}


void observer_fluxphi_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    double z_alpha, z_beta;
    double z2, phi2, eps;
    double psi_alpha_dot, psi_beta_dot, phi_dot;

    if (!s->initialized) {
        observer_init(s);
    }

    /* -----------------------------------------------------
     * 1) 用当前状态计算 z = Psi_hat - L i
     * ----------------------------------------------------- */
    z_alpha = s->psi_alpha_hat - OBS_L * i_alpha;
    z_beta  = s->psi_beta_hat  - OBS_L * i_beta;

    z2   = z_alpha * z_alpha + z_beta * z_beta;
    phi2 = s->phi_hat * s->phi_hat;

    /* 观测误差项：|z|^2 - phi_hat^2 */
    eps = z2 - phi2;

    /* -----------------------------------------------------
     * 2) 非线性磁链观测器离散更新（前向欧拉）
     * ----------------------------------------------------- */
    psi_alpha_dot = u_alpha - OBS_R * i_alpha - 2.0 * OBS_Q * z_alpha * eps;
    psi_beta_dot  = u_beta  - OBS_R * i_beta  - 2.0 * OBS_Q * z_beta  * eps;
    phi_dot       = OBS_Q * s->phi_hat * eps;

    s->psi_alpha_hat += Ts * psi_alpha_dot;
    s->psi_beta_hat  += Ts * psi_beta_dot;
    s->phi_hat       += Ts * phi_dot;

    if (s->phi_hat < OBS_PHI_MIN) {
        s->phi_hat = OBS_PHI_MIN;
    }

    /* -----------------------------------------------------
     * 3) 用更新后的状态重新计算 z
     *    后面的角度和速度，都从 z 提取
     * ----------------------------------------------------- */
    z_alpha = s->psi_alpha_hat - OBS_L * i_alpha;
    z_beta  = s->psi_beta_hat  - OBS_L * i_beta;

#if (OBS_SPEED_MODE_DEFAULT == OBS_SPEED_ATAN_DIFF)

    {
        /* =================================================
         * 方式 A：atan2 + 微分 + 一阶低通
         * 优点：简单直观
         * 缺点：速度会更敏感，低速容易抖
         * ================================================= */
        double theta_hat = atan2(z_beta, z_alpha);
        double dtheta = wrap_pm_pi(theta_hat - s->theta_prev);
        double omega_raw = dtheta / Ts;

        s->omega_e_hat_f += Ts * OBS_OMEGA_LPF_WC * (omega_raw - s->omega_e_hat_f);
        s->theta_prev = theta_hat;

        out->theta_e_hat = theta_hat;
        out->omega_e_hat = s->omega_e_hat_f;
    }

#elif (OBS_SPEED_MODE_DEFAULT == OBS_SPEED_PLL)

    {
        /* =================================================
         * 方式 B：磁链向量正交 PLL
         *
         * 构造 PLL 角度单位向量：
         *   c = cos(theta_pll)
         *   s = sin(theta_pll)
         *
         * 用 z 与其正交误差作为相位误差：
         *   e = (-z_alpha*s + z_beta*c) / |z|
         *
         * 这样：
         * - theta_pll 锁到 z 的方向
         * - omega_hat 由 PLL 环路直接输出
         * ================================================= */
        double c = cos(s->pll_theta);
        double ss = sin(s->pll_theta);
        double znorm = sqrt(z_alpha * z_alpha + z_beta * z_beta);
        double e_pll;
        double omega_hat;

        if (znorm < OBS_Z_NORM_EPS) {
            znorm = OBS_Z_NORM_EPS;
        }

        /* 正交相位误差，做归一化，减小幅值变化影响 */
        e_pll = (-z_alpha * ss + z_beta * c) / znorm;

        /* PLL PI */
        s->pll_omega_i += Ts * OBS_PLL_KI * e_pll;
        omega_hat = OBS_PLL_KP * e_pll + s->pll_omega_i;

        /* 积分出角度 */
        s->pll_theta += Ts * omega_hat;
        s->pll_theta = wrap_pm_pi(s->pll_theta);

        out->theta_e_hat = s->pll_theta;
        out->omega_e_hat = omega_hat;
    }

#else
#error "Unknown OBS_SPEED_MODE_DEFAULT"
#endif

    /* 公共输出 */
    out->psi_alpha_hat = s->psi_alpha_hat;
    out->psi_beta_hat  = s->psi_beta_hat;
    out->phi_hat       = s->phi_hat;
}