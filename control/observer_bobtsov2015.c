#include "observer_bobtsov2015.h"
#include <math.h>
#include <string.h>

/* =========================================================
 * Bobtsov et al., Automatica 2015
 *
 * 当前实现选择：
 * 1) 严格保留论文的 9 维“位置观测器”主干：式(6)~(13)
 * 2) 先只做 paper 主体，不引入 Rs 自适应
 * 3) 得到转子磁链分量 x_hat = q + eta_hat 后
 *    不再做 atan2+微分，而是直接喂给磁链向量 PLL
 * 4) 这样更适合接你现有工程里的 FOC 闭环
 *
 * 说明：
 * - 这里的 x_hat 对应论文中的 x = lambda_m * C(theta)
 * - 对 SPMSM 而言，x_hat 的方向就是电角度方向
 * - 输出给 FOC 的 theta_e_hat / omega_e_hat 全由 PLL 产生
 * ========================================================= */

#define BOB15_R             (1.025)
#define BOB15_L             (9e-3)

/* 论文实验中 α, γ 常取这个量级，先从保守值起步 */
#define BOB15_ALPHA         (100.0)
#define BOB15_GAMMA_ETA     (0.5)

/* 直接复用你主线里已经验证过的向量 PLL 风格 */
#define BOB15_PLL_KP        (300.0)
#define BOB15_PLL_KI        (30000.0)
#define BOB15_X_NORM_EPS    (1e-6)
#define BOB15_LAM_EPS       (1e-8)

static int g_init = 0;

/* ξ1..ξ9 */
static double g_xi1 = 0.0;
static double g_xi2 = 0.0;
static double g_xi3 = 0.0;
static double g_xi4 = 0.0;
static double g_xi5 = 0.0;
static double g_xi6 = 0.0;
static double g_xi7 = 0.0;
static double g_xi8 = 0.0;
static double g_xi9 = 0.0;

/* PLL 内部状态 */
static double g_pll_theta = 0.0;
static double g_pll_omega_i = 0.0;

static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < -3.14159265358979323846) x += 2.0 * 3.14159265358979323846;
    return x;
}

void observer_bobtsov2015_init(ObserverState *s)
{
    (void)s;

    g_xi1 = g_xi2 = 0.0;
    g_xi3 = g_xi4 = 0.0;
    g_xi5 = 0.0;
    g_xi6 = g_xi7 = 0.0;
    g_xi8 = g_xi9 = 0.0;

    g_pll_theta = 0.0;
    g_pll_omega_i = 0.0;

    g_init = 1;
}

void observer_bobtsov2015_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}

void observer_bobtsov2015_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    double q1, q2;
    double q_norm2;
    double Omega1, Omega2;
    double yreg;
    double err_eta;
    double xhat_alpha, xhat_beta;
    double lam_hat_alpha, lam_hat_beta;
    double lam_mag;
    double c, ss, xnorm, e_pll, omega_hat;

    (void)s;

    if (!g_init) {
        observer_bobtsov2015_init(NULL);
    }

    /* ---------------------------------------------
     * 1) 论文式(6): ξ12dot = v, ξ34dot = i
     * --------------------------------------------- */
    g_xi1 += Ts * u_alpha;
    g_xi2 += Ts * u_beta;
    g_xi3 += Ts * i_alpha;
    g_xi4 += Ts * i_beta;

    /* ---------------------------------------------
     * 2) 论文式(7): q = ξ12 - R ξ34 - L i
     *    q 实际上就是 rotor-flux 分量 x 的“去偏前”版本
     * --------------------------------------------- */
    q1 = g_xi1 - BOB15_R * g_xi3 - BOB15_L * i_alpha;
    q2 = g_xi2 - BOB15_R * g_xi4 - BOB15_L * i_beta;
    q_norm2 = q1 * q1 + q2 * q2;

    /* ---------------------------------------------
     * 3) 论文式(8)~(11): 形成 y, Omega
     *    ξ5dot  = -α(ξ5 + |q|^2)
     *    ξ67dot = -α(ξ67 - 2q)
     *    y      = -α(|q|^2 + ξ5)
     *    Ω      = α(2q - ξ67)
     * --------------------------------------------- */
    g_xi5 += Ts * (-BOB15_ALPHA * (g_xi5 + q_norm2));
    g_xi6 += Ts * (-BOB15_ALPHA * (g_xi6 - 2.0 * q1));
    g_xi7 += Ts * (-BOB15_ALPHA * (g_xi7 - 2.0 * q2));

    yreg   = -BOB15_ALPHA * (q_norm2 + g_xi5);
    Omega1 =  BOB15_ALPHA * (2.0 * q1 - g_xi6);
    Omega2 =  BOB15_ALPHA * (2.0 * q2 - g_xi7);

    /* ---------------------------------------------
     * 4) 论文式(12): η 梯度估计
     *    ξ89dot = Γ Ω (y - Ω^T ξ89)
     *    这里 Γ 先取 γI
     * --------------------------------------------- */
    err_eta = yreg - (Omega1 * g_xi8 + Omega2 * g_xi9);
    g_xi8  += Ts * (BOB15_GAMMA_ETA * Omega1 * err_eta);
    g_xi9  += Ts * (BOB15_GAMMA_ETA * Omega2 * err_eta);

    /* ---------------------------------------------
     * 5) x_hat = q + eta_hat
     *    这是转子磁链分量估计，方向就是电角度方向
     * --------------------------------------------- */
    xhat_alpha = q1 + g_xi8;
    xhat_beta  = q2 + g_xi9;

    /* 总磁链估计 lambda_hat = L i + x_hat */
    lam_hat_alpha = BOB15_L * i_alpha + xhat_alpha;
    lam_hat_beta  = BOB15_L * i_beta  + xhat_beta;
    lam_mag = sqrt(xhat_alpha * xhat_alpha + xhat_beta * xhat_beta);
    if (lam_mag < BOB15_LAM_EPS) lam_mag = BOB15_LAM_EPS;

    /* ---------------------------------------------
     * 6) 用转子磁链向量 x_hat 进 PLL
     *    这一步就是你建议的做法
     * --------------------------------------------- */
    c = cos(g_pll_theta);
    ss = sin(g_pll_theta);
    xnorm = sqrt(xhat_alpha * xhat_alpha + xhat_beta * xhat_beta);
    if (xnorm < BOB15_X_NORM_EPS) {
        xnorm = BOB15_X_NORM_EPS;
    }

    e_pll = (-xhat_alpha * ss + xhat_beta * c) / xnorm;

    g_pll_omega_i += Ts * BOB15_PLL_KI * e_pll;
    omega_hat = BOB15_PLL_KP * e_pll + g_pll_omega_i;

    g_pll_theta += Ts * omega_hat;
    g_pll_theta = wrap_pm_pi(g_pll_theta);

    /* ---------------------------------------------
     * 7) 输出
     * --------------------------------------------- */
    memset(out, 0, sizeof(*out));

    out->theta_e_hat = g_pll_theta;
    out->omega_e_hat = omega_hat;

    /* 为了沿用现有 logger，这里仍把总磁链写到 psi_hat */
    out->psi_alpha_hat = lam_hat_alpha;
    out->psi_beta_hat  = lam_hat_beta;
    out->phi_hat       = lam_mag;

    out->eta1_hat = g_xi8;
    out->eta2_hat = g_xi9;

    /* 借现有调试槽位把 paper 关键中间量打出来 */
    out->q1 = q1;
    out->q2 = q2;
    out->yreg = yreg;
    out->detQ = Omega1 * Omega1 + Omega2 * Omega2;
    out->xi1 = xhat_alpha;
    out->xi2 = xhat_beta;
    out->z21 = Omega1;
    out->z22 = Omega2;
}