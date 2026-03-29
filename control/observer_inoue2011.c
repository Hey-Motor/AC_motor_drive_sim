#include "observer_inoue2011.h"
#include "../config/config.h"
#include <math.h>
#include <string.h>

/* =========================================================
 * Inoue 2011 / Morimoto 2002 observer
 *
 * 一、无感本体
 * ---------------------------------------------------------
 * - Morimoto 2002：
 *   在估计旋转 gamma-delta 坐标系中，用 extended-EMF 模型建立
 *   least-order observer；先做交叉耦合补偿，再分别估 e_gamma/e_delta。
 *
 *   论文关系式（对应代码中的实现口径）：
 *   (A) v_gamma1 = v_gamma + omega_hat * Lq * i_delta
 *   (B) v_delta1 = v_delta - omega_hat * Ld * i_gamma
 *   (C) dot(z_gamma) = -g*z_gamma + g*v_gamma1 + (g^2*Ld - g*R)*i_gamma
 *       e_gamma_hat  = z_gamma - g*Ld*i_gamma
 *   (D) dot(z_delta) = -g*z_delta + g*v_delta1 + (g^2*Lq - g*R)*i_delta
 *       e_delta_hat  = z_delta - g*Lq*i_delta
 *   (E) theta_err_hat = atan2(-e_gamma_hat, e_delta_hat)    [Scheme A]
 *   (F) omega_hat_raw = K_ep * theta_err_hat + xi
 *       dot(xi) = K_ei * theta_err_hat
 *       dot(theta_hat) = omega_hat_raw
 *
 * 二、在线 Rs 辨识
 * ---------------------------------------------------------
 * - Inoue 2011：使用 q-axis identification model；在无感系统里，
 *   用 estimated gamma-delta frame 的变量替代 dq 变量。
 * - 这里采用标量 RLS：
 *     Y = v_delta(k-1) - Ld*omega(k-1)*i_gamma(k-1) - psi_f*omega(k-1)
 *     Z = i_delta(k-1)
 *     Theta = R
 *
 * 三、工程取舍
 * ---------------------------------------------------------
 * - 当前接口未把 Params 传进 observer，因此这里先把电机参数写成宏，
 *   与你工程当前默认参数保持一致。后续若你要做正式版，可再改成从 params
 *   或 observer 专用结构体中取值。
 * - 当前分支优先保证：
 *   1) 与现有 observer_if / observer_select 接口兼容；
 *   2) 先能单独跑 observer；
 *   3) Rs 辨识可通过宏单独开关。
 * ========================================================= */

/* =========================
 * 电机参数（按当前工程默认值）
 * 若你切回 IPMSM，请直接改这里
 * ========================= */
#if (USE_IPMSM == 1)
#define INOUE_LD_INIT            (9.0e-3)
#define INOUE_LQ_INIT            (18.50e-3)   /* 这里给个占位，按你的实际电机再改 */
#define INOUE_PSI_F_INIT         (0.249)
#else
#define INOUE_LD_INIT            (9.0e-3)
#define INOUE_LQ_INIT            (9.0e-3)
#define INOUE_PSI_F_INIT         (0.249)
#endif

#define INOUE_R_INIT             (0.3)

/* =========================
 * Morimoto 2002 观察器 / 速度修正起始参数
 * 论文实验：g=600 rad/s, wn=45 rad/s, zeta=0.5, omega LPF=100 rad/s
 * 这些是“论文电机”的起始量级，不保证你当前电机无需再调。
 * ========================= */
#define INOUE_GAMMA_OBS          (600.0)
#define INOUE_PLL_WN             (100.0)
#define INOUE_PLL_ZETA           (0.5)
#define INOUE_OMEGA_LPF_WC       (100.0)

#define INOUE_K_EP               (2.0 * INOUE_PLL_ZETA * INOUE_PLL_WN)
#define INOUE_K_EI               (INOUE_PLL_WN * INOUE_PLL_WN)

/* =========================
 * Rs 在线辨识开关与参数
 * ========================= */
#define INOUE_ENABLE_RS_RLS      (1)
#define INOUE_RLS_LAMBDA         (0.9995)
#define INOUE_RLS_P0             (1000.0)
#define INOUE_RLS_MIN_EXC_I      (0.20)      /* |i_delta| 太小时，不更新 */
#define INOUE_RLS_MIN_OMEGA      (20.0)      /* 低于该电角速度时，不更新 */
#define INOUE_RLS_DI_MAX         (10.0)      /* delta 轴电流变化过快时，不更新（A/s） */
#define INOUE_R_MIN              (0.05)
#define INOUE_R_MAX              (20.0)

/* =========================
 * 数值保护
 * ========================= */
#define INOUE_EPS                (1e-12)
#define INOUE_TWO_PI             (6.28318530717958647692)
#define INOUE_PI                 (3.14159265358979323846)

static double wrap_pm_pi(double x)
{
    while (x > INOUE_PI)  x -= INOUE_TWO_PI;
    while (x < -INOUE_PI) x += INOUE_TWO_PI;
    return x;
}

static double clampd(double x, double xmin, double xmax)
{
    if (x < xmin) return xmin;
    if (x > xmax) return xmax;
    return x;
}

/* =========================
 * file-static 内部状态
 * ========================= */
static double g_theta_hat;
static double g_omega_hat_raw;
static double g_omega_hat_f;
static double g_xi_pll;

static double g_z_gamma;
static double g_z_delta;
static double g_e_gamma_hat;
static double g_e_delta_hat;

static double g_R_hat;
static double g_P_rls;

/* 上一拍 gamma-delta / RLS 所需缓存 */
static int    g_has_prev;
static double g_prev_i_gamma;
static double g_prev_i_delta;
static double g_prev_v_gamma;
static double g_prev_v_delta;
static double g_prev_omega_hat_f;

/* 仅供外部调试接口兼容，当前不使用 */
void observer_inoue2011_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}

void observer_inoue2011_init(ObserverState *s)
{
    if (s != NULL) {
        memset(s, 0, sizeof(*s));
        s->initialized = 1;
    }

    g_theta_hat      = 0.0;
    g_omega_hat_raw  = 0.0;
    g_omega_hat_f    = 0.0;
    g_xi_pll         = 0.0;

    g_z_gamma        = 0.0;
    g_z_delta        = 0.0;
    g_e_gamma_hat    = 0.0;
    g_e_delta_hat    = 0.0;

    g_R_hat          = INOUE_R_INIT;
    g_P_rls          = INOUE_RLS_P0;

    g_has_prev       = 0;
    g_prev_i_gamma   = 0.0;
    g_prev_i_delta   = 0.0;
    g_prev_v_gamma   = 0.0;
    g_prev_v_delta   = 0.0;
    g_prev_omega_hat_f = 0.0;
}

void observer_inoue2011_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    const double Ld = INOUE_LD_INIT;
    const double Lq = INOUE_LQ_INIT;
    const double psi_f = INOUE_PSI_F_INIT;
    const double g = INOUE_GAMMA_OBS;
    const double R_use = (INOUE_ENABLE_RS_RLS ? g_R_hat : INOUE_R_INIT);

    double c, sn;
    double i_gamma, i_delta;
    double v_gamma, v_delta;
    double v_gamma1, v_delta1;
    double dz_gamma, dz_delta;
    double theta_err_hat;
    double e_alpha_hat, e_beta_hat;

    if ((s == NULL) || (out == NULL)) {
        return;
    }
    if (!s->initialized) {
        observer_inoue2011_init(s);
    }

    /* -----------------------------------------------------
     * 1) alpha-beta -> estimated gamma-delta
     * ----------------------------------------------------- */
    c  = cos(g_theta_hat);
    sn = sin(g_theta_hat);

    i_gamma =  c * i_alpha + sn * i_beta;
    i_delta = -sn * i_alpha + c * i_beta;

    v_gamma =  c * u_alpha + sn * u_beta;
    v_delta = -sn * u_alpha + c * u_beta;

    /* -----------------------------------------------------
     * 2) Morimoto 2002: 交叉耦合补偿
     * ----------------------------------------------------- */
    v_gamma1 = v_gamma + g_omega_hat_f * Lq * i_delta;
    v_delta1 = v_delta - g_omega_hat_f * Ld * i_gamma;

    /* -----------------------------------------------------
     * 3) least-order observer for extended EMF
     * ----------------------------------------------------- */
    dz_gamma = -g * g_z_gamma
             +  g * v_gamma1
             + (g * g * Ld - g * R_use) * i_gamma;
    g_z_gamma += Ts * dz_gamma;
    g_e_gamma_hat = g_z_gamma - g * Ld * i_gamma;

    dz_delta = -g * g_z_delta
             +  g * v_delta1
             + (g * g * Lq - g * R_use) * i_delta;
    g_z_delta += Ts * dz_delta;
    g_e_delta_hat = g_z_delta - g * Lq * i_delta;

    /* -----------------------------------------------------
     * 4) Scheme A: 由 extended EMF 求位置误差
     * ----------------------------------------------------- */
    theta_err_hat = atan2(-g_e_gamma_hat, g_e_delta_hat + INOUE_EPS);
    theta_err_hat = wrap_pm_pi(theta_err_hat);

    /* -----------------------------------------------------
     * 5) 位置 / 速度 PI 修正
     * ----------------------------------------------------- */
    g_xi_pll += Ts * INOUE_K_EI * theta_err_hat;
    g_omega_hat_raw = INOUE_K_EP * theta_err_hat + g_xi_pll;

    g_omega_hat_f += Ts * INOUE_OMEGA_LPF_WC * (g_omega_hat_raw - g_omega_hat_f);

    g_theta_hat += Ts * g_omega_hat_raw;
    g_theta_hat = wrap_pm_pi(g_theta_hat);

#if (INOUE_ENABLE_RS_RLS == 1)
    /* -----------------------------------------------------
     * 6) Inoue 2011: q-axis identification model for Rs
     *    注意：这里用 estimated gamma-delta frame 的变量替代 dq 变量。
     *    采用上一拍量，贴近论文的离散形式。
     * ----------------------------------------------------- */
    if (g_has_prev) {
        double Y, Z, den, K, err;
        double di_delta_dt;
        int excite_ok;

        di_delta_dt = 0;//(i_delta - g_prev_i_delta) / (Ts + INOUE_EPS);
        excite_ok = (fabs(g_prev_i_delta) >= INOUE_RLS_MIN_EXC_I)
                 && (fabs(g_prev_omega_hat_f) >= INOUE_RLS_MIN_OMEGA)
                 && (fabs(di_delta_dt) <= INOUE_RLS_DI_MAX);

        if (excite_ok) {
            Y = g_prev_v_delta
              - Ld * g_prev_omega_hat_f * g_prev_i_gamma
              - psi_f * g_prev_omega_hat_f;
            Z = g_prev_i_delta;

            den = INOUE_RLS_LAMBDA + Z * g_P_rls * Z;
            K   = (g_P_rls * Z) / (den + INOUE_EPS);
            err = Y - Z * g_R_hat;

            g_R_hat += K * err;
            g_R_hat  = clampd(g_R_hat, INOUE_R_MIN, INOUE_R_MAX);

            g_P_rls = (1.0 / INOUE_RLS_LAMBDA) * (g_P_rls - K * Z * g_P_rls);
            if (g_P_rls < INOUE_EPS) {
                g_P_rls = INOUE_EPS;
            }

            out->q1 = Y;
            out->q2 = Z;
            out->q3 = err;
            out->q4 = K;
            out->q5 = den;
            out->q6 = di_delta_dt;
        }
    }
#endif

    /* 更新上一拍缓存 */
    g_prev_i_gamma     = i_gamma;
    g_prev_i_delta     = i_delta;
    g_prev_v_gamma     = v_gamma;
    g_prev_v_delta     = v_delta;
    g_prev_omega_hat_f = g_omega_hat_f;
    g_has_prev         = 1;

    /* -----------------------------------------------------
     * 7) 调试输出：
     *    - theta/omega：给 FOC 用
     *    - psi_alpha/beta_hat：这里复用为 alpha-beta 下的 extended EMF
     *    - eta1/eta2：复用为 gamma-delta 下的 extended EMF
     *    - phi_hat：位置误差
     *    - beta_hat：未低通的速度
     *    - detQ：复用为 RLS 协方差 P
     * ----------------------------------------------------- */
    e_alpha_hat =  c * g_e_gamma_hat - sn * g_e_delta_hat;
    e_beta_hat  =  sn * g_e_gamma_hat + c * g_e_delta_hat;

    out->theta_e_hat = g_theta_hat;
    out->omega_e_hat = g_omega_hat_f;

    out->psi_alpha_hat = e_alpha_hat;
    out->psi_beta_hat  = e_beta_hat;
    out->phi_hat       = theta_err_hat;

    out->R_hat         = g_R_hat;
    out->eta1_hat      = g_e_gamma_hat;
    out->eta2_hat      = g_e_delta_hat;
    out->beta_hat      = g_omega_hat_raw;
    out->detQ          = g_P_rls;

    /* 这些槽位在本分支里没有特定物理含义，清零即可 */
    out->yreg = 0.0;
    out->z21  = 0.0;
    out->z22  = 0.0;
    out->xi1  = g_z_gamma;
    out->xi2  = g_z_delta;

    /* 同步写回公共 state，方便你后续 logger / 调试统一看 */
    s->psi_alpha_hat = e_alpha_hat;
    s->psi_beta_hat  = e_beta_hat;
    s->phi_hat       = theta_err_hat;
    s->theta_prev    = g_theta_hat;
    s->omega_e_hat_f = g_omega_hat_f;
    s->pll_theta     = g_theta_hat;
    s->pll_omega_i   = g_xi_pll;
}
