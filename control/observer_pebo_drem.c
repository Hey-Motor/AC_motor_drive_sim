/* =========================================================
 * EXPERIMENTAL BRANCH
 * 当前版本仅用于 open-loop observer-only 调试
 * 不建议直接接管 FOC 闭环
 * 明天再基于“只估 eta1 eta2”做简化版
 * ========================================================= */

#include "observer_pebo_drem.h"
#include <math.h>
#include <string.h>

/* =========================================================
 * Bazylev, Pyrkin, Bobtsov 2018
 * "Position and speed observer for PMSM with unknown stator resistance"
 *
 * 论文核心结构：
 * 1) z1dot = v, z2dot = i
 * 2) 构造 6 维回归 -> DREM -> 估计 R, eta1, eta2
 * 3) lambda_hat = z1 - R_hat z2 + eta_hat
 * 4) x_hat = lambda_hat - L i
 * 5) omega_hat = psi1 + beta_hat
 * 6) theta_hat = psi2 + theta0_hat
 * 7) theta0_hat 用混合动态反演，而不是 atan
 *
 * 当前实现说明：
 * - DREM 的 5 个附加算子：使用 5 个 delay operator
 * - 位置直接输出 electrical angle（给 FOC 更自然）
 * - 内部状态采用 file-static，改动最小；当前工程只有一个 observer 实例
 * ========================================================= */

/* =========================
 * 已知电机参数（对应论文假设 A4）
 * 请与你当前 SPMSM plant 参数保持一致
 * ========================= */
#define PEBO_L              (9e-3)
#define PEBO_J              (0.002379)
#define PEBO_B              (1e-4)
#define PEBO_NP             (3.0)

/* =========================
 * DREM / observer 参数
 * ========================= */
#define PEBO_ALPHA          (100.0)     /* αp/(p+α) */
#define PEBO_GAMMA_R        (100.0)
#define PEBO_GAMMA_ETA1     (100.0)
#define PEBO_GAMMA_ETA2     (100.0)

#define PEBO_ZETA           (500.0)     /* β观测器滤波参数 ς */
#define PEBO_GAMMA_BETA     (400.0)

#define PEBO_GAMMA_K1       (50.0)      /* theta0 混合恢复增益 */
#define PEBO_GAMMA_K2       (50.0)

/* 5 个 delay operator，对应 DREM 的 Hk */
/* DREM 的 5 个 LTI 滤波器参数
 * Hk(p) = eps_k / (p + eps_k)
 * 按论文仿真参数设置
 */
#define PEBO_EPS1           (100.0)
#define PEBO_EPS2           (180.0)
#define PEBO_EPS3           (260.0)
#define PEBO_EPS4           (340.0)
#define PEBO_EPS5           (420.0)
#define PEBO_DET_EPS        (1e-10)
#define PEBO_LAM_EPS        (1e-8)

/* 每个 LTI 滤波器都要对 6 个 q 和 1 个 y 做一阶滤波 */
static double g_qf[5][6];
static double g_yf[5];

/* =========================
 * file-static 内部状态
 * 这样不需要大改 ObserverState
 * 若以后你要多 observer 实例，再搬入结构体
 * ========================= */
static int g_init = 0;

/* z1dot=v, z2dot=i */
static double g_z1[2];
static double g_z2[2];

/* α/(p+α) 的低通内部状态，用于得到 αp/(p+α) */
static double g_lp_raw[7];


/* 参数估计 */
static double g_R_hat = 0.1;
static double g_eta1_hat = 0.0;
static double g_eta2_hat = 0.0;

/* 磁链 / x */
static double g_lambda_hat[2];
static double g_x_hat[2];
static double g_lambda_m_hat = 0.0;

/* ψ1, β, ψ2 */
static double g_psi1_mech = 0.0;
static double g_beta_hat = 0.0;
static double g_psi2_elec = 0.0;

/* β 观测器滤波状态 */
static double g_lp_i[2];
static double g_lp_w[2];
static double g_lp_psi3[2];

/* θ0 混合恢复内部状态（直接用 electrical angle 形式） */
static double g_k1e = 0.1;
static double g_k2e = 0.1;


/* =========================================================
 * 基础工具
 * ========================================================= */
static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < -3.14159265358979323846) x += 2.0 * 3.14159265358979323846;
    return x;
}

static double wrap_0_2pi(double x)
{
    while (x >= 2.0 * 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < 0.0) x += 2.0 * 3.14159265358979323846;
    return x;
}

static double sign1(double x)
{
    return (x >= 0.0) ? 1.0 : -1.0;
}

/* 一阶低通: a/(p+a) */
static double lpf_step(double *x, double a, double u, double Ts)
{
    *x += Ts * a * (u - *x);
    return *x;
}

/* 一阶 LTI 滤波器
 * H(p) = eps / (p + eps)
 */
static double drem_filter_step(double *x, double eps, double u, double Ts)
{
    *x += Ts * eps * (u - *x);
    return *x;
}

/* αp/(p+α) = α*(u - α/(p+α)u 的输出) */
static double hp_like_step(double *xlp, double a, double u, double Ts)
{
    double ylp = lpf_step(xlp, a, u, Ts);
    return a * (u - ylp);
}


/* 6x6 高斯消元：同时给 det 和解 x，供 DREM 使用 */
static int solve6(double A[6][6], double b[6], double x[6], double *det_out)
{
    int i, j, k, piv;
    double det = 1.0;
    double aug[6][7];

    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) aug[i][j] = A[i][j];
        aug[i][6] = b[i];
    }

    for (i = 0; i < 6; ++i) {
        double maxv = fabs(aug[i][i]);
        piv = i;
        for (k = i + 1; k < 6; ++k) {
            double v = fabs(aug[k][i]);
            if (v > maxv) { maxv = v; piv = k; }
        }

        if (maxv < PEBO_DET_EPS) {
            *det_out = 0.0;
            return 0;
        }

        if (piv != i) {
            for (j = i; j < 7; ++j) {
                double tmp = aug[i][j];
                aug[i][j] = aug[piv][j];
                aug[piv][j] = tmp;
            }
            det = -det;
        }

        det *= aug[i][i];

        for (k = i + 1; k < 6; ++k) {
            double f = aug[k][i] / aug[i][i];
            for (j = i; j < 7; ++j) {
                aug[k][j] -= f * aug[i][j];
            }
        }
    }

    for (i = 5; i >= 0; --i) {
        double s = aug[i][6];
        for (j = i + 1; j < 6; ++j) s -= aug[i][j] * x[j];
        x[i] = s / aug[i][i];
    }

    *det_out = det;
    return 1;
}

/* electrical θ0 的区间投影 */
static double project_k1(double k1e, double xbar2)
{
    k1e = wrap_0_2pi(k1e);

    /* A1=(0,pi), A2=(pi,2pi) */
    if (xbar2 >= 0.0) {
        if (k1e <= 0.0) k1e = 1e-6;
        if (k1e >= 3.14159265358979323846) k1e = 3.14159265358979323846 - 1e-6;
    } else {
        if (k1e <= 3.14159265358979323846) k1e = 3.14159265358979323846 + 1e-6;
        if (k1e >= 2.0 * 3.14159265358979323846) k1e = 2.0 * 3.14159265358979323846 - 1e-6;
    }
    return k1e;
}

static double project_k2(double k2e, double xbar1)
{
    const double pi = 3.14159265358979323846;
    const double hpi = 0.5 * pi;
    const double thpi = 1.5 * pi;

    k2e = wrap_0_2pi(k2e);

    /* A3=[0,pi/2) U (3pi/2,2pi), A4=(pi/2,3pi/2) */
    if (xbar1 >= 0.0) {
        if (k2e > hpi && k2e < thpi) {
            if (k2e < pi) k2e = hpi - 1e-6;
            else          k2e = thpi + 1e-6;
        }
    } else {
        if (k2e <= hpi) k2e = hpi + 1e-6;
        if (k2e >= thpi) k2e = thpi - 1e-6;
    }

    return wrap_0_2pi(k2e);
}


/* =========================================================
 * init / debug hook
 * ========================================================= */
void observer_pebo_drem_init(ObserverState *s)
{
    (void)s;

    memset(g_z1, 0, sizeof(g_z1));
    memset(g_z2, 0, sizeof(g_z2));
    memset(g_lp_raw, 0, sizeof(g_lp_raw));
    memset(g_lp_i, 0, sizeof(g_lp_i));
    memset(g_lp_w, 0, sizeof(g_lp_w));
    memset(g_lp_psi3, 0, sizeof(g_lp_psi3));
    memset(g_qf, 0, sizeof(g_qf));
    memset(g_yf, 0, sizeof(g_yf));

    g_R_hat = 0.1;
    g_eta1_hat = 0.0;
    g_eta2_hat = 0.0;

    g_lambda_hat[0] = 0.0;
    g_lambda_hat[1] = 0.0;
    g_x_hat[0] = 0.0;
    g_x_hat[1] = 0.0;
    g_lambda_m_hat = 0.0;

    g_psi1_mech = 0.0;
    g_beta_hat = 0.0;
    g_psi2_elec = 0.0;

    g_k1e = 0.1;
    g_k2e = 0.1;

    g_init = 1;
}

void observer_pebo_drem_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}


/* =========================================================
 * 主更新
 * ========================================================= */
void observer_pebo_drem_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    double xi1, xi2;
    double raw[6], q[6], yreg;
    double Qe[6][6], Ye[6], mu_inst[6], detQ = 0.0;
    double tau_e, omega_hat_mech, omega_hat_elec;
    double jx0, jx1, di0, di1, w0, w1, wf0, wf1, psi30, psi31, yhat0, yhat1, ebeta;
    double xbar1, xbar2, c2, s2, theta0e_hat;
    int ready = 0;

    (void)s;

    if (!g_init) observer_pebo_drem_init(NULL);

    /* -----------------------------------------------------
     * 1) z1dot = v, z2dot = i
     * ----------------------------------------------------- */
    g_z1[0] += Ts * u_alpha;
    g_z1[1] += Ts * u_beta;
    g_z2[0] += Ts * i_alpha;
    g_z2[1] += Ts * i_beta;

    /* xi = z1 - L i */
    xi1 = g_z1[0] - PEBO_L * i_alpha;
    xi2 = g_z1[1] - PEBO_L * i_beta;

    /* -----------------------------------------------------
     * 2) 构造原始 6 维回归
     * 来自论文式(12)
     * 为了工程离散实现，这里统一通过 αp/(p+α) 形成
     * ----------------------------------------------------- */
    raw[0] = -2.0 * (g_z2[0] * xi1 + g_z2[1] * xi2);
    raw[1] =  2.0 * xi1;
    raw[2] =  2.0 * xi2;
    raw[3] =  (g_z2[0] * g_z2[0] + g_z2[1] * g_z2[1]);
    raw[4] = -2.0 * g_z2[0];
    raw[5] = -2.0 * g_z2[1];

    q[0] = hp_like_step(&g_lp_raw[0], PEBO_ALPHA, raw[0], Ts);
    q[1] = hp_like_step(&g_lp_raw[1], PEBO_ALPHA, raw[1], Ts);
    q[2] = hp_like_step(&g_lp_raw[2], PEBO_ALPHA, raw[2], Ts);
    q[3] = hp_like_step(&g_lp_raw[3], PEBO_ALPHA, raw[3], Ts);
    q[4] = hp_like_step(&g_lp_raw[4], PEBO_ALPHA, raw[4], Ts);
    q[5] = hp_like_step(&g_lp_raw[5], PEBO_ALPHA, raw[5], Ts);

    yreg = hp_like_step(&g_lp_raw[6], PEBO_ALPHA, -(xi1 * xi1 + xi2 * xi2), Ts);

   /* =====================================================
 * DREM: 原始回归 + 5 个 LTI 滤波回归
 * Hk(p) = eps_k / (p + eps_k)
 * ===================================================== */
{
    int j;
    const double eps_list[5] = {
        PEBO_EPS1, PEBO_EPS2, PEBO_EPS3, PEBO_EPS4, PEBO_EPS5
    };

    /* 第 0 行：原始回归 */
    for (j = 0; j < 6; ++j) {
        Qe[0][j] = q[j];
    }
    Ye[0] = yreg;

    /* 第 1~5 行：5 个 LTI 滤波后的回归 */
    for (int k = 0; k < 5; ++k) {
        for (j = 0; j < 6; ++j) {
            g_qf[k][j] = drem_filter_step(&g_qf[k][j], eps_list[k], q[j], Ts);
            Qe[k + 1][j] = g_qf[k][j];
        }

        g_yf[k] = drem_filter_step(&g_yf[k], eps_list[k], yreg, Ts);
        Ye[k + 1] = g_yf[k];
    }

    ready = 1;
}

    if (ready) {
        double A[6][6];
        double b[6];
        int i, j;

        for (i = 0; i < 6; ++i) {
            for (j = 0; j < 6; ++j) A[i][j] = Qe[i][j];
            b[i] = Ye[i];
        }

        /* 用 det(Qe) 和即时解 mu_inst 形成 Yl = phi * mu_l */
        if (solve6(A, b, mu_inst, &detQ)) {
            double phi = detQ;
            double Y1 = phi * mu_inst[0];
            double Y2 = phi * mu_inst[1];
            double Y3 = phi * mu_inst[2];

            g_R_hat    += Ts * PEBO_GAMMA_R    * phi * (Y1 - phi * g_R_hat);
            g_eta1_hat += Ts * PEBO_GAMMA_ETA1 * phi * (Y2 - phi * g_eta1_hat);
            g_eta2_hat += Ts * PEBO_GAMMA_ETA2 * phi * (Y3 - phi * g_eta2_hat);
        }
    }

    /* -----------------------------------------------------
     * 4) 磁链观测器：
     * lambda_hat = z1 - R_hat z2 + eta_hat
     * ----------------------------------------------------- */
    g_lambda_hat[0] = g_z1[0] - g_R_hat * g_z2[0] + g_eta1_hat;
    g_lambda_hat[1] = g_z1[1] - g_R_hat * g_z2[1] + g_eta2_hat;

    /* x_hat = lambda_hat - L i */
    g_x_hat[0] = g_lambda_hat[0] - PEBO_L * i_alpha;
    g_x_hat[1] = g_lambda_hat[1] - PEBO_L * i_beta;

    g_lambda_m_hat = sqrt(g_x_hat[0] * g_x_hat[0] + g_x_hat[1] * g_x_hat[1]);
    if (g_lambda_m_hat < PEBO_LAM_EPS) g_lambda_m_hat = PEBO_LAM_EPS;

    /* -----------------------------------------------------
     * 5) 速度观测器：
     * tau_e = np * i^T J lambda
     * psi1 = 1/(Jp+B) tau_e
     * beta_hat 用论文式(34)-(36)
     * ----------------------------------------------------- */
    tau_e = -PEBO_NP * (i_alpha * g_lambda_hat[1] - i_beta * g_lambda_hat[0]);

    /* ψ1 机械角速度分量 */
    g_psi1_mech += Ts * ((tau_e - PEBO_B * g_psi1_mech) / PEBO_J);

    /* np J x_hat */
    jx0 = -PEBO_NP * g_x_hat[1];
    jx1 =  PEBO_NP * g_x_hat[0];

    /* ψ3 = ς/(p+ς)(np J x_hat) */
    psi30 = lpf_step(&g_lp_psi3[0], PEBO_ZETA, jx0, Ts);
    psi31 = lpf_step(&g_lp_psi3[1], PEBO_ZETA, jx1, Ts);

    /* L * ςp/(p+ς) i */
    {
        double low_i0 = lpf_step(&g_lp_i[0], PEBO_ZETA, i_alpha, Ts);
        double low_i1 = lpf_step(&g_lp_i[1], PEBO_ZETA, i_beta, Ts);
        di0 = PEBO_ZETA * (i_alpha - low_i0);
        di1 = PEBO_ZETA * (i_beta  - low_i1);
    }

    /* ς/(p+ς)(-Rhat i + v - ψ1 npJxhat) */
    w0 = -g_R_hat * i_alpha + u_alpha - g_psi1_mech * jx0;
    w1 = -g_R_hat * i_beta  + u_beta  - g_psi1_mech * jx1;

    wf0 = lpf_step(&g_lp_w[0], PEBO_ZETA, w0, Ts);
    wf1 = lpf_step(&g_lp_w[1], PEBO_ZETA, w1, Ts);

    yhat0 = PEBO_L * di0 - wf0;
    yhat1 = PEBO_L * di1 - wf1;

    ebeta = psi30 * (yhat0 - g_beta_hat * psi30)
          + psi31 * (yhat1 - g_beta_hat * psi31);

    g_beta_hat += Ts * PEBO_GAMMA_BETA * ebeta;

    omega_hat_mech = g_psi1_mech + g_beta_hat;
    omega_hat_elec = PEBO_NP * omega_hat_mech;

    /* -----------------------------------------------------
     * 6) 位置观测器：
     * theta_hat = psi2 + theta0_hat
     * 这里直接用 electrical angle 形式
     * ----------------------------------------------------- */
    g_psi2_elec += Ts * omega_hat_elec;
    g_psi2_elec = wrap_0_2pi(g_psi2_elec);

    /* xbar = R(-psi2_elec) x_hat / lambda_m_hat */
    c2 = cos(g_psi2_elec);
    s2 = sin(g_psi2_elec);
    xbar1 = ( c2 * g_x_hat[0] + s2 * g_x_hat[1]) / g_lambda_m_hat;
    xbar2 = (-s2 * g_x_hat[0] + c2 * g_x_hat[1]) / g_lambda_m_hat;

    /* 论文式(47)-(49) 的 electrical-angle 等价实现 */
    g_k1e += Ts * PEBO_GAMMA_K1 * (-xbar1 + cos(g_k1e)) * sign1(xbar2);
    g_k1e  = project_k1(g_k1e, xbar2);

    g_k2e += Ts * PEBO_GAMMA_K2 * ( xbar2 - sin(g_k2e)) * sign1(xbar1);
    g_k2e  = project_k2(g_k2e, xbar1);

    theta0e_hat = (fabs(xbar1) <= fabs(xbar2)) ? g_k1e : g_k2e;

    // out->theta_e_hat = wrap_0_2pi(g_psi2_elec + theta0e_hat);
    // out->omega_e_hat = omega_hat_elec;

    /* 调试用：先绕开混合 theta0 恢复，只看 x_hat 的方向是否正确 */
    out->theta_e_hat = atan2(g_x_hat[1], g_x_hat[0]);

    /* 速度仍然保留当前 beta/psi1 链 */
    out->omega_e_hat = omega_hat_elec;

    /* 这几个保留给 logger 用 */
    out->psi_alpha_hat = g_lambda_hat[0];
    out->psi_beta_hat  = g_lambda_hat[1];
    out->phi_hat       = g_lambda_m_hat;

    out->R_hat    = g_R_hat;
    out->eta1_hat = g_eta1_hat;
    out->eta2_hat = g_eta2_hat;
    out->beta_hat = g_beta_hat;
    out->detQ     = detQ;
    out->q1 = q[0];
    out->q2 = q[1];
    out->q3 = q[2];
    out->q4 = q[3];
    out->q5 = q[4];
    out->q6 = q[5];
    out->yreg = yreg;
    
    out->z21 = g_z2[0];
    out->z22 = g_z2[1];
    out->xi1 = xi1;
    out->xi2 = xi2;
    /* 如果你后面想在 logger 里直接看这些量，
     * 建议把 ObserverOutput 扩展出对应字段：
     *   R_hat, eta1_hat, eta2_hat, beta_hat, drem_phi
     * 当前先不强制改 types.h，避免侵入过大
     */
}