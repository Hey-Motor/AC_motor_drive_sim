#include "observer_bernard2021.h"
#include "../config/config.h"
#include <math.h>
#include <string.h>

/* =========================================================
 * Bernard & Praly 2021: first practical integration version
 *
 * 一、论文核心对象
 * ---------------------------------------------------------
 * 论文在 alpha-beta 固定坐标下，以 total flux Ψ 为状态，
 * 采用：
 *   Ψdot = u - R i
 *   |Ψ - L i|^2 = Φ^2
 *
 * 然后对三个不同 λ 构造 Luenberger 变换 T_λ，得到 21 个滤波器状态：
 *   a_λ, b_λ(2), c_λ(2), d_λ, e_λ
 * 并通过静态方程：
 *   M(R,t) χ(R,t) = rhs(R,t)
 *   J(R,t) = 0
 * 来恢复 Ψ 与 R。
 *
 * 二、当前工程里的实现取舍
 * ---------------------------------------------------------
 * 1) 动态部分：严格对应论文式(26a)~(26e)；
 * 2) 静态部分：为了便于并入你现有工程，先做 fixed-grid 的 Argmin |J| 搜索，
 *    并叠加 iq 符号筛选。它非常接近论文 IV-A，但没有做“精确根检测 + 多解全返回”；
 * 3) 速度：实现论文式(38)，而不是简单 atan 微分；
 * 4) 启动初期滤波器还没忘记初值时，用一个轻量的电压模型磁链积分兜底，
 *    避免一上电 theta_hat 完全空缺。
 *
 * 三、务必说明
 * ---------------------------------------------------------
 * - 这是“可运行工程版”，不是把论文的所有定理细节一字不差搬进 C。
 * - 尤其静态部分，这里优先选了：工程可用 / 接口兼容 / 便于日志调试。
 * - 如果后面你要做“论文完全对口版”，再继续把：
 *   1) root sign-change 精细检测
 *   2) 多候选同时输出
 *   3) saliency 修正 J_tilde / Theta_s
 *   这些分支继续加进去。
 * ========================================================= */

/* =========================
 * 有效电机参数（论文主体按 SPMSM）
 * ========================= */
#if (USE_IPMSM == 1)
/* 若你仍在 IPMSM 工况下强行试验，本实现只能看成近似版。 */
#define BER21_L_INIT              (18.5e-3)
#else
#define BER21_L_INIT              (9.0e-3)
#endif
#define BER21_PSI_F_INIT          (0.249)
#define BER21_R_INIT              (1.025)

/* =========================
 * 三组 lambda：论文理想数据示例用 20/30/40
 * ========================= */
#define BER21_LAMBDA1             (20.0)
#define BER21_LAMBDA2             (30.0)
#define BER21_LAMBDA3             (40.0)

/* =========================
 * 静态部分：固定网格搜索参数
 * ========================= */
#define BER21_R_MIN               (0.05)
#define BER21_R_MAX               (5.00)
#define BER21_R_STEP              (0.01)
#define BER21_R_UPDATE_DT         (0.005)   /* 5 ms 更新一次 R */
#define BER21_WAIT_TIME           (0.50)    /* 论文示例：0.5 s 后开始估计 */
#define BER21_R_INIT_LOCK_TIME    (0.30)    /* 启动后先在名义初值附近锁分支 */
#define BER21_R_INIT_LOCK_G       (0.25)
#define BER21_R_LOCAL_G           (0.20)    /* 平时采用小幅移动网格 */
#define BER21_R_RESCUE_PERIOD     (1e9)     /* 当前版本禁用周期救援 */
#define BER21_J_BAD_RATIO         (1e9)     /* 当前版本禁用“J 恶化触发救援” */
#define BER21_R_TRACK_ALPHA       (0.08)    /* 平滑输出 Rs，只用于日志/观测 */
#define BER21_DETM_EPS            (1e-10)
#define BER21_IQ_SIGN_TARGET      (+1.0)    /* 默认按电动工况筛选 iq>0 */
#define BER21_IQ_SIGN_EPS         (0.05)

/* =========================
 * 速度估计器：论文式(38)
 * ideal data 示例给的是 rho=1000, k=500
 * ========================= */
#define BER21_SPEED_RHO           (1000.0)
#define BER21_SPEED_K             (2000.0)
#define BER21_OMEGA_MAX           (4000.0)

/* =========================
 * 启动兜底：轻量电压模型磁链积分
 * ========================= */
#define BER21_VM_LEAK_WC          (5.0)

#define BER21_PI                  (3.14159265358979323846)
#define BER21_TWO_PI              (6.28318530717958647692)
#define BER21_EPS                 (1e-12)

#define BER21_OMEGA_WC            (2.0 * BER21_PI * 20.0)  /* 速度低通截止频率，先取20Hz */
#define BER21_OMEGA_RAW_MAX       (2000.0)                 /* 微分后原始速度限幅 */
#define BER21_USE_ALT_SPEED       (1)                      /* 1=用替代速度法 */

typedef struct {
    double a;
    double b_alpha;
    double b_beta;
    double c_alpha;
    double c_beta;
    double d;
    double e;
    double lambda;
} Ber21LamState;

static Ber21LamState g_f[3];

static double g_time;
static double g_R_hat;
static double g_R_branch;
static double g_theta_hat;
static double g_omega_hat;

/* 论文式(38) 的内部状态 */
static double g_chi_s_hat;
static double g_chi_c_hat;
static double g_nu_hat;

/* 启动兜底用的简单电压模型总磁链 */
static double g_vm_psi_alpha;
static double g_vm_psi_beta;

/* 调试输出缓存 */
static double g_last_detM;
static double g_last_Jabs;
static double g_last_rotor_flux_alpha;
static double g_last_rotor_flux_beta;
static double g_last_Psi_alpha;
static double g_last_Psi_beta;
static double g_last_best_alt_R;
static double g_last_best_alt_J;
static double g_last_best_iq;
static double g_last_x3a;
static double g_last_x3b;

static int g_initialized;

static double g_theta_prev_speed = 0.0;
static double g_omega_lpf = 0.0;
static int    g_speed_inited = 0;

static double wrap_pm_pi(double x)
{
    while (x > BER21_PI)  x -= BER21_TWO_PI;
    while (x < -BER21_PI) x += BER21_TWO_PI;
    return x;
}

static double clampd(double x, double xmin, double xmax)
{
    if (x < xmin) return xmin;
    if (x > xmax) return xmax;
    return x;
}

static double sign_target_ok(double iq)
{
    if (fabs(iq) < BER21_IQ_SIGN_EPS) {
        return 1.0;
    }
    return (iq * BER21_IQ_SIGN_TARGET >= 0.0) ? 1.0 : 0.0;
}

static void ber21_filter_init(Ber21LamState *s, double lambda)
{
    memset(s, 0, sizeof(*s));
    s->lambda = lambda;
}

static void ber21_filter_step(
    Ber21LamState *s,
    double Ts,
    double L,
    double Phi,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta)
{
    double lam = s->lambda;
    double b_dot_alpha, b_dot_beta;
    double c_dot_alpha, c_dot_beta;
    double a_dot, d_dot, e_dot;
    double ci, bu, bi, cu;
    double i_norm2;

    ci = s->c_alpha * i_alpha + s->c_beta * i_beta;
    bu = s->b_alpha * u_alpha + s->b_beta * u_beta;
    bi = s->b_alpha * i_alpha + s->b_beta * i_beta;
    cu = s->c_alpha * u_alpha + s->c_beta * u_beta;
    i_norm2 = i_alpha * i_alpha + i_beta * i_beta;

    a_dot = -lam * (s->a - ci + bu);

    b_dot_alpha = -lam * (s->b_alpha - 2.0 * i_alpha);
    b_dot_beta  = -lam * (s->b_beta  - 2.0 * i_beta);

    c_dot_alpha = -lam * (s->c_alpha + 2.0 * u_alpha + 2.0 * lam * L * i_alpha);
    c_dot_beta  = -lam * (s->c_beta  + 2.0 * u_beta  + 2.0 * lam * L * i_beta);

    d_dot = -lam * (s->d - bi);

    e_dot = -lam * (s->e - cu + lam * lam * L * L * i_norm2 - lam * lam * Phi * Phi);

    s->a       += Ts * a_dot;
    s->b_alpha += Ts * b_dot_alpha;
    s->b_beta  += Ts * b_dot_beta;
    s->c_alpha += Ts * c_dot_alpha;
    s->c_beta  += Ts * c_dot_beta;
    s->d       += Ts * d_dot;
    s->e       += Ts * e_dot;
}

static int ber21_eval_candidate(
    double Rcand,
    double i_alpha,
    double i_beta,
    double *Psi_alpha,
    double *Psi_beta,
    double *rotor_flux_alpha,
    double *rotor_flux_beta,
    double *theta,
    double *iq_cand,
    double *detM,
    double *Jabs,
    double *x3a,
    double *x3b)
{
    const double l1 = g_f[0].lambda;
    const double l2 = g_f[1].lambda;
    const double l3 = g_f[2].lambda;

    double q1, q2, q3;
    double C1a, C1b, C2a, C2b, C3a, C3b;
    double M11, M12, M21, M22, det;
    double rhs1, rhs2;
    double chi_a, chi_b;
    double T1, T2, T3;
    double J;
    double rf_a, rf_b;
    double th, iqv;

    /* qk = e_k - a_k R - d_k R^2 */
    q1 = g_f[0].e - g_f[0].a * Rcand - g_f[0].d * Rcand * Rcand;
    q2 = g_f[1].e - g_f[1].a * Rcand - g_f[1].d * Rcand * Rcand;
    q3 = g_f[2].e - g_f[2].a * Rcand - g_f[2].d * Rcand * Rcand;

    /* Ck = c_k + R b_k */
    C1a = g_f[0].c_alpha + Rcand * g_f[0].b_alpha;
    C1b = g_f[0].c_beta  + Rcand * g_f[0].b_beta;
    C2a = g_f[1].c_alpha + Rcand * g_f[1].b_alpha;
    C2b = g_f[1].c_beta  + Rcand * g_f[1].b_beta;
    C3a = g_f[2].c_alpha + Rcand * g_f[2].b_alpha;
    C3b = g_f[2].c_beta  + Rcand * g_f[2].b_beta;

    /*
     * M_lambda = [[ l2^2, -l1^2, 0 ],
     *             [ 0,    l3^2, -l2^2 ]]
     * M = M_lambda * Lambda * [C1;C2;C3]
     */
    M11 = (l2 * l2) * (l1 * C1a) - (l1 * l1) * (l2 * C2a);
    M12 = (l2 * l2) * (l1 * C1b) - (l1 * l1) * (l2 * C2b);
    M21 = (l3 * l3) * (l2 * C2a) - (l2 * l2) * (l3 * C3a);
    M22 = (l3 * l3) * (l2 * C2b) - (l2 * l2) * (l3 * C3b);

    det = M11 * M22 - M12 * M21;
    if (detM != NULL) *detM = det;

    /* det(M)=A R^2 + B R + C，供日志观察“坏点”位置 */
    if ((x3a != NULL) || (x3b != NULL)) {
        double b1a = g_f[0].b_alpha, b1b = g_f[0].b_beta;
        double b2a = g_f[1].b_alpha, b2b = g_f[1].b_beta;
        double b3a = g_f[2].b_alpha, b3b = g_f[2].b_beta;
        double c1a = g_f[0].c_alpha, c1b = g_f[0].c_beta;
        double c2a = g_f[1].c_alpha, c2b = g_f[1].c_beta;
        double c3a = g_f[2].c_alpha, c3b = g_f[2].c_beta;
        double A11 = (l2 * l2) * (l1 * b1a) - (l1 * l1) * (l2 * b2a);
        double A12 = (l2 * l2) * (l1 * b1b) - (l1 * l1) * (l2 * b2b);
        double A21 = (l3 * l3) * (l2 * b2a) - (l2 * l2) * (l3 * b3a);
        double A22 = (l3 * l3) * (l2 * b2b) - (l2 * l2) * (l3 * b3b);
        double B11 = (l2 * l2) * (l1 * c1a) - (l1 * l1) * (l2 * c2a);
        double B12 = (l2 * l2) * (l1 * c1b) - (l1 * l1) * (l2 * c2b);
        double B21 = (l3 * l3) * (l2 * c2a) - (l2 * l2) * (l3 * c3a);
        double B22 = (l3 * l3) * (l2 * c2b) - (l2 * l2) * (l3 * c3b);
        double Aq = A11 * A22 - A12 * A21;
        double Bq = A11 * B22 + B11 * A22 - A12 * B21 - B12 * A21;
        double Cq = B11 * B22 - B12 * B21;
        double disc = Bq * Bq - 4.0 * Aq * Cq;
        if (fabs(Aq) < BER21_EPS) {
            if (fabs(Bq) < BER21_EPS) {
                if (x3a) *x3a = 0.0;
                if (x3b) *x3b = 0.0;
            } else {
                if (x3a) *x3a = -Cq / Bq;
                if (x3b) *x3b = -Cq / Bq;
            }
        } else if (disc >= 0.0) {
            double sq = sqrt(disc);
            if (x3a) *x3a = (-Bq - sq) / (2.0 * Aq);
            if (x3b) *x3b = (-Bq + sq) / (2.0 * Aq);
        } else {
            if (x3a) *x3a = 0.0;
            if (x3b) *x3b = 0.0;
        }
    }

    if (fabs(det) < BER21_DETM_EPS) {
        return 0;
    }

    rhs1 = (l2 * l2) * q1 - (l1 * l1) * q2;
    rhs2 = (l3 * l3) * q2 - (l2 * l2) * q3;

    chi_a = ( rhs1 * M22 - M12 * rhs2) / det;
    chi_b = (-rhs1 * M21 + M11 * rhs2) / det;

    T1 = l1 * l1 * (chi_a * chi_a + chi_b * chi_b)
       + l1 * ((C1a * chi_a) + (C1b * chi_b))
       + g_f[0].a * Rcand + g_f[0].d * Rcand * Rcand - g_f[0].e;

    T2 = l2 * l2 * (chi_a * chi_a + chi_b * chi_b)
       + l2 * ((C2a * chi_a) + (C2b * chi_b))
       + g_f[1].a * Rcand + g_f[1].d * Rcand * Rcand - g_f[1].e;

    T3 = l3 * l3 * (chi_a * chi_a + chi_b * chi_b)
       + l3 * ((C3a * chi_a) + (C3b * chi_b))
       + g_f[2].a * Rcand + g_f[2].d * Rcand * Rcand - g_f[2].e;

    J = l1 * l1 * (g_f[0].e - T1)
      + l2 * l2 * (g_f[1].e - T2)
      + l3 * l3 * (g_f[2].e - T3);

    rf_a = chi_a - BER21_L_INIT * i_alpha;
    rf_b = chi_b - BER21_L_INIT * i_beta;
    th = atan2(rf_b, rf_a + BER21_EPS);
    iqv = -sin(th) * i_alpha + cos(th) * i_beta;

    if (Psi_alpha) *Psi_alpha = chi_a;
    if (Psi_beta)  *Psi_beta  = chi_b;
    if (rotor_flux_alpha) *rotor_flux_alpha = rf_a;
    if (rotor_flux_beta)  *rotor_flux_beta  = rf_b;
    if (theta) *theta = th;
    if (iq_cand) *iq_cand = iqv;
    if (Jabs) *Jabs = fabs(J);

    return 1;
}

void observer_bernard2021_set_debug_truth(double theta_e_true, double omega_e_true)
{
    (void)theta_e_true;
    (void)omega_e_true;
}

void observer_bernard2021_init(ObserverState *s)
{
    if (s != NULL) {
        memset(s, 0, sizeof(*s));
        s->initialized = 1;
    }

    ber21_filter_init(&g_f[0], BER21_LAMBDA1);
    ber21_filter_init(&g_f[1], BER21_LAMBDA2);
    ber21_filter_init(&g_f[2], BER21_LAMBDA3);

    g_time = 0.0;
    g_R_hat = BER21_R_INIT;
    g_R_branch = BER21_R_INIT;
    g_theta_hat = 0.0;
    g_omega_hat = 0.0;

    g_chi_s_hat = 0.0;
    g_chi_c_hat = 0.0;
    g_nu_hat = 0.0;

    g_vm_psi_alpha = 0.0;
    g_vm_psi_beta  = 0.0;

    g_last_detM = 0.0;
    g_last_Jabs = 0.0;
    g_last_rotor_flux_alpha = 0.0;
    g_last_rotor_flux_beta  = 0.0;
    g_last_Psi_alpha = 0.0;
    g_last_Psi_beta  = 0.0;
    g_last_best_alt_R = 0.0;
    g_last_best_alt_J = 0.0;
    g_last_best_iq = 0.0;
    g_last_x3a = 0.0;
    g_last_x3b = 0.0;

    g_initialized = 1;
}

void observer_bernard2021_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out)
{
    double Psi_alpha = 0.0;
    double Psi_beta  = 0.0;
    double rf_alpha = 0.0;
    double rf_beta  = 0.0;
    double theta_flux;
    double iq_here;
    double rho = BER21_SPEED_RHO;
    double k = BER21_SPEED_K;
    double theta_for_speed;
    double omega_speed;

    if ((s == NULL) || (out == NULL)) {
        return;
    }
    if (!g_initialized || !s->initialized) {
        observer_bernard2021_init(s);
    }

    memset(out, 0, sizeof(*out));

    /* -----------------------------------------------------
     * 1) 论文式(26)：动态滤波器部分
     * ----------------------------------------------------- */
    ber21_filter_step(&g_f[0], Ts, BER21_L_INIT, BER21_PSI_F_INIT,
                      u_alpha, u_beta, i_alpha, i_beta);
    ber21_filter_step(&g_f[1], Ts, BER21_L_INIT, BER21_PSI_F_INIT,
                      u_alpha, u_beta, i_alpha, i_beta);
    ber21_filter_step(&g_f[2], Ts, BER21_L_INIT, BER21_PSI_F_INIT,
                      u_alpha, u_beta, i_alpha, i_beta);

    /* -----------------------------------------------------
     * 2) 简单电压模型兜底：用于滤波器尚未进入稳态阶段
     * ----------------------------------------------------- */
    g_vm_psi_alpha += Ts * (u_alpha - g_R_branch * i_alpha - BER21_VM_LEAK_WC * g_vm_psi_alpha);
    g_vm_psi_beta  += Ts * (u_beta  - g_R_branch * i_beta  - BER21_VM_LEAK_WC * g_vm_psi_beta);

    /* -----------------------------------------------------
     * 3) 静态部分：每隔 dt_R 在固定网格上搜索最优 R
     *    这里先采用“Argmin |J| + iq 符号优先”的工程版。
     * ----------------------------------------------------- */
    g_time += Ts;
    if (g_time >= BER21_WAIT_TIME) {
        static double r_timer = 0.0;
        static double rescue_timer = 0.0;
        r_timer += Ts;
        rescue_timer += Ts;

        if (r_timer >= BER21_R_UPDATE_DT) {
            double best_R_sign = g_R_branch;
            double best_J_sign = 1e300;
            double best_R_any  = g_R_branch;
            double best_J_any  = 1e300;
            double best_det_sign = 0.0;
            double best_det_any  = 0.0;
            double best_theta_sign = g_theta_hat;
            double best_theta_any  = g_theta_hat;
            double best_Psi_a_sign = 0.0, best_Psi_b_sign = 0.0;
            double best_rf_a_sign = 0.0, best_rf_b_sign = 0.0;
            double best_Psi_a_any = 0.0,  best_Psi_b_any = 0.0;
            double best_rf_a_any = 0.0,   best_rf_b_any = 0.0;
            double best_iq_sign = 0.0,    best_iq_any = 0.0;
            double best_x3a = 0.0, best_x3b = 0.0;
            double r;
            double r_lo, r_hi;
            int used_global_rescue = 0;

            /* 论文的工程口径：小幅网格围绕初值/当前估计移动；
             * 只有在局部判据明显变坏时，才做一次全局救援搜索。 */
            if (g_time < BER21_WAIT_TIME + BER21_R_INIT_LOCK_TIME) {
                r_lo = BER21_R_INIT - BER21_R_INIT_LOCK_G;
                r_hi = BER21_R_INIT + BER21_R_INIT_LOCK_G;
            } else {
                r_lo = g_R_branch - BER21_R_LOCAL_G;
                r_hi = g_R_branch + BER21_R_LOCAL_G;
            }
            if (r_lo < BER21_R_MIN) r_lo = BER21_R_MIN;
            if (r_hi > BER21_R_MAX) r_hi = BER21_R_MAX;

            for (r = r_lo; r <= r_hi + 0.5 * BER21_R_STEP; r += BER21_R_STEP) {
                double psi_a_c, psi_b_c, rf_a_c, rf_b_c, theta_c, iq_c, detM_c, Jabs_c, xa, xb;
                int ok;

                ok = ber21_eval_candidate(r, i_alpha, i_beta,
                                          &psi_a_c, &psi_b_c,
                                          &rf_a_c, &rf_b_c,
                                          &theta_c, &iq_c,
                                          &detM_c, &Jabs_c,
                                          &xa, &xb);
                if (!ok) {
                    continue;
                }

                if (Jabs_c < best_J_any) {
                    best_J_any = Jabs_c;
                    best_R_any = r;
                    best_det_any = detM_c;
                    best_theta_any = theta_c;
                    best_Psi_a_any = psi_a_c;
                    best_Psi_b_any = psi_b_c;
                    best_rf_a_any = rf_a_c;
                    best_rf_b_any = rf_b_c;
                    best_iq_any = iq_c;
                    best_x3a = xa;
                    best_x3b = xb;
                }

                if (sign_target_ok(iq_c) > 0.5 && Jabs_c < best_J_sign) {
                    best_J_sign = Jabs_c;
                    best_R_sign = r;
                    best_det_sign = detM_c;
                    best_theta_sign = theta_c;
                    best_Psi_a_sign = psi_a_c;
                    best_Psi_b_sign = psi_b_c;
                    best_rf_a_sign = rf_a_c;
                    best_rf_b_sign = rf_b_c;
                    best_iq_sign = iq_c;
                }
            }

            {
                double Jref = (g_last_Jabs > 1e-12) ? g_last_Jabs : 1e-12;
                int near_boundary = 0;
                int local_bad = (best_J_any > BER21_J_BAD_RATIO * Jref);
                int periodic_rescue = (rescue_timer >= BER21_R_RESCUE_PERIOD);

                if (near_boundary || local_bad || periodic_rescue) {
                    used_global_rescue = 1;
                    best_R_sign = g_R_branch;
                    best_J_sign = 1e300;
                    best_R_any  = g_R_branch;
                    best_J_any  = 1e300;
                    best_det_sign = 0.0;
                    best_det_any  = 0.0;
                    best_theta_sign = g_theta_hat;
                    best_theta_any  = g_theta_hat;
                    best_Psi_a_sign = 0.0; best_Psi_b_sign = 0.0;
                    best_rf_a_sign = 0.0;  best_rf_b_sign = 0.0;
                    best_Psi_a_any = 0.0;  best_Psi_b_any = 0.0;
                    best_rf_a_any = 0.0;   best_rf_b_any = 0.0;
                    best_iq_sign = 0.0;    best_iq_any = 0.0;
                    best_x3a = 0.0; best_x3b = 0.0;

                    for (r = BER21_R_MIN; r <= BER21_R_MAX + 0.5 * BER21_R_STEP; r += BER21_R_STEP) {
                        double psi_a_c, psi_b_c, rf_a_c, rf_b_c, theta_c, iq_c, detM_c, Jabs_c, xa, xb;
                        int ok;

                        ok = ber21_eval_candidate(r, i_alpha, i_beta,
                                                  &psi_a_c, &psi_b_c,
                                                  &rf_a_c, &rf_b_c,
                                                  &theta_c, &iq_c,
                                                  &detM_c, &Jabs_c,
                                                  &xa, &xb);
                        if (!ok) {
                            continue;
                        }

                        if (Jabs_c < best_J_any) {
                            best_J_any = Jabs_c;
                            best_R_any = r;
                            best_det_any = detM_c;
                            best_theta_any = theta_c;
                            best_Psi_a_any = psi_a_c;
                            best_Psi_b_any = psi_b_c;
                            best_rf_a_any = rf_a_c;
                            best_rf_b_any = rf_b_c;
                            best_iq_any = iq_c;
                            best_x3a = xa;
                            best_x3b = xb;
                        }

                        if (sign_target_ok(iq_c) > 0.5 && Jabs_c < best_J_sign) {
                            best_J_sign = Jabs_c;
                            best_R_sign = r;
                            best_det_sign = detM_c;
                            best_theta_sign = theta_c;
                            best_Psi_a_sign = psi_a_c;
                            best_Psi_b_sign = psi_b_c;
                            best_rf_a_sign = rf_a_c;
                            best_rf_b_sign = rf_b_c;
                            best_iq_sign = iq_c;
                        }
                    }
                    rescue_timer = 0.0;
                }
            }

            {
                double R_sel;
                double theta_sel;
                double det_sel;
                double J_sel;
                double Psi_a_sel, Psi_b_sel;
                double rf_a_sel, rf_b_sel;
                double iq_sel;

                if (best_J_sign < 1e299) {
                    R_sel     = best_R_sign;
                    theta_sel = best_theta_sign;
                    det_sel   = best_det_sign;
                    J_sel     = best_J_sign;
                    Psi_a_sel = best_Psi_a_sign;
                    Psi_b_sel = best_Psi_b_sign;
                    rf_a_sel  = best_rf_a_sign;
                    rf_b_sel  = best_rf_b_sign;
                    iq_sel    = best_iq_sign;
                } else if (best_J_any < 1e299) {
                    R_sel     = best_R_any;
                    theta_sel = best_theta_any;
                    det_sel   = best_det_any;
                    J_sel     = best_J_any;
                    Psi_a_sel = best_Psi_a_any;
                    Psi_b_sel = best_Psi_b_any;
                    rf_a_sel  = best_rf_a_any;
                    rf_b_sel  = best_rf_b_any;
                    iq_sel    = best_iq_any;
                } else {
                    R_sel     = g_R_branch;
                    theta_sel = g_theta_hat;
                    det_sel   = g_last_detM;
                    J_sel     = g_last_Jabs;
                    Psi_a_sel = g_last_Psi_alpha;
                    Psi_b_sel = g_last_Psi_beta;
                    rf_a_sel  = g_last_rotor_flux_alpha;
                    rf_b_sel  = g_last_rotor_flux_beta;
                    iq_sel    = g_last_best_iq;
                }

                /* 原始分支值：给下一次搜索/分支连续性使用 */
                g_R_branch = R_sel;

                /* 平滑输出值：仅用于日志与外部观测，避免段内台阶跳变 */
                g_R_hat += BER21_R_TRACK_ALPHA * (g_R_branch - g_R_hat);

                g_theta_hat = theta_sel;
                g_last_detM = det_sel;
                g_last_Jabs = J_sel;
                g_last_Psi_alpha = Psi_a_sel;
                g_last_Psi_beta  = Psi_b_sel;
                g_last_rotor_flux_alpha = rf_a_sel;
                g_last_rotor_flux_beta  = rf_b_sel;
                g_last_best_iq = iq_sel;
            }

            g_last_best_alt_R = best_R_any;
            g_last_best_alt_J = best_J_any;
            g_last_x3a = best_x3a;
            g_last_x3b = best_x3b;
            (void)used_global_rescue;

            r_timer = 0.0;
        }

        /* 用当前原始分支值连续计算 Ψ_hat / theta_hat */
        if (ber21_eval_candidate(g_R_branch, i_alpha, i_beta,
                                 &Psi_alpha, &Psi_beta,
                                 &rf_alpha, &rf_beta,
                                 &theta_flux, &iq_here,
                                 &g_last_detM, &g_last_Jabs,
                                 &g_last_x3a, &g_last_x3b)) {
            g_theta_hat = theta_flux;
            g_last_Psi_alpha = Psi_alpha;
            g_last_Psi_beta  = Psi_beta;
            g_last_rotor_flux_alpha = rf_alpha;
            g_last_rotor_flux_beta  = rf_beta;
            g_last_best_iq = iq_here;
        } else {
            rf_alpha = g_vm_psi_alpha - BER21_L_INIT * i_alpha;
            rf_beta  = g_vm_psi_beta  - BER21_L_INIT * i_beta;
            g_theta_hat = atan2(rf_beta, rf_alpha + BER21_EPS);
            g_last_Psi_alpha = g_vm_psi_alpha;
            g_last_Psi_beta  = g_vm_psi_beta;
            g_last_rotor_flux_alpha = rf_alpha;
            g_last_rotor_flux_beta  = rf_beta;
            g_last_best_iq = -sin(g_theta_hat) * i_alpha + cos(g_theta_hat) * i_beta;
        }
    } else {
        rf_alpha = g_vm_psi_alpha - BER21_L_INIT * i_alpha;
        rf_beta  = g_vm_psi_beta  - BER21_L_INIT * i_beta;
        g_theta_hat = atan2(rf_beta, rf_alpha + BER21_EPS);
        g_last_Psi_alpha = g_vm_psi_alpha;
        g_last_Psi_beta  = g_vm_psi_beta;
        g_last_rotor_flux_alpha = rf_alpha;
        g_last_rotor_flux_beta  = rf_beta;
        g_last_best_iq = -sin(g_theta_hat) * i_alpha + cos(g_theta_hat) * i_beta;
        g_last_detM = 0.0;
        g_last_Jabs = 0.0;
    }

    g_theta_hat = wrap_pm_pi(g_theta_hat);

    /* -----------------------------------------------------
     * 4) 论文式(38)：速度估计器
     * ----------------------------------------------------- */
    // theta_for_speed = g_theta_hat;
    // g_chi_s_hat += Ts * (-rho * g_chi_s_hat - (rho * rho + g_nu_hat) * sin(theta_for_speed));
    // g_chi_c_hat += Ts * (-rho * g_chi_c_hat - (rho * rho + g_nu_hat) * cos(theta_for_speed));
    // g_nu_hat    += Ts * (k * ((g_chi_c_hat * cos(theta_for_speed) + g_chi_s_hat * sin(theta_for_speed)) + rho));

    // omega_speed = g_chi_s_hat * cos(theta_for_speed) - g_chi_c_hat * sin(theta_for_speed);
    // g_omega_hat = clampd(omega_speed, -BER21_OMEGA_MAX, BER21_OMEGA_MAX);
    /* -----------------------------------------------------
    * 4) 速度估计
    *    theta_hat 仍来自 Bernard 的转子磁链方向
    *    omega_hat 暂时改为：角度解缠微分 + 一阶低通
    * ----------------------------------------------------- */
    {
        double theta_for_speed = g_theta_hat;

    #if (BER21_USE_ALT_SPEED == 1)

        if (!g_speed_inited) {
            g_theta_prev_speed = theta_for_speed;
            g_omega_lpf = 0.0;
            g_omega_hat = 0.0;
            g_speed_inited = 1;
        } else {
            double dtheta = wrap_pm_pi(theta_for_speed - g_theta_prev_speed);
            double omega_raw = dtheta / Ts;
            double a = BER21_OMEGA_WC * Ts;

            if (a > 1.0) a = 1.0;

            if (omega_raw >  BER21_OMEGA_RAW_MAX) omega_raw =  BER21_OMEGA_RAW_MAX;
            if (omega_raw < -BER21_OMEGA_RAW_MAX) omega_raw = -BER21_OMEGA_RAW_MAX;

            g_omega_lpf += a * (omega_raw - g_omega_lpf);
            g_omega_hat = g_omega_lpf;

            g_theta_prev_speed = theta_for_speed;
        }

        /* 这三个量先保留一个合理值，避免日志或别处引用时乱掉 */
        g_chi_s_hat = 0.0;
        g_chi_c_hat = 0.0;
        g_nu_hat    = g_omega_hat;

    #else

        /* 原论文式(38) */
        g_chi_s_hat += Ts * (-rho * g_chi_s_hat - (rho * rho + g_nu_hat) * sin(theta_for_speed));
        g_chi_c_hat += Ts * (-rho * g_chi_c_hat - (rho * rho + g_nu_hat) * cos(theta_for_speed));
        g_nu_hat    += Ts * (k * ((g_chi_c_hat * cos(theta_for_speed) + g_chi_s_hat * sin(theta_for_speed)) + rho));

        omega_speed = g_chi_s_hat * cos(theta_for_speed) - g_chi_c_hat * sin(theta_for_speed);
        g_omega_hat = omega_speed;

    #endif
    }

    /* -----------------------------------------------------
     * 5) 输出映射：尽量贴合你现有 logger 的槽位
     * ----------------------------------------------------- */
    out->theta_e_hat = g_theta_hat;
    out->omega_e_hat = g_omega_hat;

    /* 总磁链 Ψ_hat */
    out->psi_alpha_hat = g_last_Psi_alpha;
    out->psi_beta_hat  = g_last_Psi_beta;

    /* rotor flux x_hat = Ψ_hat - L i；幅值理想上接近 psi_f */
    out->phi_hat  = sqrt(g_last_rotor_flux_alpha * g_last_rotor_flux_alpha
                       + g_last_rotor_flux_beta  * g_last_rotor_flux_beta);

    out->R_hat    = g_R_hat;
    out->eta1_hat = g_last_rotor_flux_alpha;
    out->eta2_hat = g_last_rotor_flux_beta;
    out->beta_hat = g_last_Jabs;
    out->detQ     = g_last_detM;

    /* 调试量：
     * q1,q2 : 当前选中候选对应的 rotor flux alpha-beta
     * q3,q4 : 速度估计器 chi_s/chi_c
     * q5    : nu_hat
     * q6    : 选中候选的 iq
     * yreg  : 次优/全局最小 |J|
     * z21,z22 : det(M)=0 的两个“坏点”近似位置
     * xi1,xi2 : 当前全局最优 R 与当前选中 R
     */
    out->q1  = g_last_rotor_flux_alpha;
    out->q2  = g_last_rotor_flux_beta;
    out->q3  = g_chi_s_hat;
    out->q4  = g_chi_c_hat;
    out->q5  = g_nu_hat;
    out->q6  = g_last_best_iq;
    out->yreg = g_last_best_alt_J;

    out->z21 = g_last_x3a;
    out->z22 = g_last_x3b;
    out->xi1 = g_R_hat;      /* 平滑后的输出值 */
    out->xi2 = g_R_branch;   /* 原始分支值 */

    s->psi_alpha_hat = out->psi_alpha_hat;
    s->psi_beta_hat  = out->psi_beta_hat;
    s->phi_hat       = out->phi_hat;
    s->theta_prev    = g_theta_hat;
    s->omega_e_hat_f = g_omega_hat;
    s->pll_theta     = g_theta_hat;
    s->pll_omega_i   = g_nu_hat;
}
