#include "observer_ipmsm_nlrs.h"
#include <math.h>
#include <string.h>
#include "foc_controller.h"
/* =========================================================
 * Lin / Liu / Liang 2023
 * Joint Estimation of Resistance and Speed With Auxiliary State
 * of Active Flux for Sensorless IPMSM Drives
 *
 * 这一版实现的是论文中的“非线性参数化全阶自适应观测器”主干：
 *   - 新坐标状态 x = [psi_sigma; chi]
 *   - 参数 theta = [omega_e, Rs]
 *   - 观测器/参数律/滤波回归器 对应式(14)
 *
 * 工程说明：
 * 1) 这是一个“最小侵入、可跑”的离散化版本，采用前向欧拉。
 * 2) 当前工程 observer 接口只有 u,i 输入，因此本文件内部用 file-static 状态。
 * 3) 位置提取不直接照抄 PDF 中排版有点乱的式(18)字符形式，
 *    而是等价地先恢复 active EMF 向量 e_active_hat，再旋转得到
 *    active-flux 方向向量，最后 atan2 得到电角度。
 * 4) 这篇方法是针对 IPMSM 的，要求 Lq 已知且较准确；Ld 不参与估计。
 * ========================================================= */

/* ===== 与当前工程 IPMSM 参数保持一致 ===== */
#define OBS_LQ              (18.5e-3)
#define OBS_RS_INIT         (1.025)//(1.025)

/* ===== 论文 Appendix B 推荐量级 =====
 * gamma = diag(1.5e3, 1e4), lambda1 = 1e3, lambda2 = 2e4
 * 先按这个量级起步。
 * ===== */
#define OBS_GAMMA_W         (9.0e9) //9.0e9    2026 03 23 改之前速度估计好用的参数 rad/s 100
#define OBS_GAMMA_R         (7.0e6) //7.0e6    2026 03 23 改之前速度估计好用的参数 rad/s 100
#define OBS_LAMBDA1         (1.0e3) //1.0e3    2026 03 23 改之前速度估计好用的参数 rad/s 100
#define OBS_LAMBDA2         (2.0e4) //2.0e4    2026 03 23 改之前速度估计好用的参数 rad/s 100

/* 速度接近零时，位置公式里需要一个可靠符号 */
#define OBS_SIGN_W_THR      (5.0)
#define OBS_OMEGA_INIT      (10.0)
#define OBS_RS_MIN          (0.05)
#define OBS_RS_MAX          (10.0)
#define OBS_OMEGA_MAX       (4000.0)
#define OBS_EPS             (1e-9)

#define OBS_OMEGA_INIT_E   (0.0)

#define OBS_NU_R_SCALE   (1)
#define OBS_R_ADAPT_LPF_W   (100.0)   /* rad/s，先用这个 */
#define OBS_NU_R_UP_SCALE   (1.0)
static double g_R_adapt_raw_f = 0.0;
static int g_init = 0;

/* x_hat = [psi_sigma_hat(2); chi_hat(2)] */
static double g_psis_hat[2];
static double g_chi_hat[2];

/* theta_hat = [omega_e_hat, Rs_hat] */
static double g_omega_hat_e = OBS_OMEGA_INIT;
static double g_R_hat = OBS_RS_INIT;
static double g_omega_hat_e_obs = 0;

/* ν ∈ R^{4x2}; 列 0 对应 omega, 列 1 对应 Rs */
static double g_nu[4][2];

/* 低速时保留上一次可靠符号，避免 sign(omega_hat)=0 */
static double g_sign_omega = 1.0;

static double g_dbg_theta_true_e = 0.0;
static double g_dbg_omega_true_e = 0.0;

static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < -3.14159265358979323846) x += 2.0 * 3.14159265358979323846;
    return x;
}

static double clamp(double x, double xmin, double xmax)
{
    if (x > xmax) return xmax;
    if (x < xmin) return xmin;
    return x;
}

/* J * [a b]^T = [-b a]^T */
static void J_mul(double a, double b, double *ja, double *jb)
{
    *ja = -b;
    *jb =  a;
}

void observer_ipmsm_nlrs_init(ObserverState *s)
{
    (void)s;

    memset(g_psis_hat, 0, sizeof(g_psis_hat));
    memset(g_chi_hat, 0, sizeof(g_chi_hat));
    memset(g_nu, 0, sizeof(g_nu));

    //g_omega_hat_e = OBS_OMEGA_INIT;
    g_omega_hat_e = OBS_OMEGA_INIT_E;
    g_R_hat = OBS_RS_INIT;
    g_sign_omega = 1.0;
    g_omega_hat_e_obs  = 0;

    g_init = 1;
    g_R_adapt_raw_f = 0.0;
}

void observer_ipmsm_nlrs_set_debug_truth(double theta_e_true, double omega_e_true)
{
    g_dbg_theta_true_e = theta_e_true;
    g_dbg_omega_true_e = omega_e_true;
}

void observer_ipmsm_nlrs_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    double y_psis_a, y_psis_b;
    double eps_a, eps_b;
    double Ji_a, Ji_b;
    double JuRi_a, JuRi_b;
    double dphi_domega[4];
    double dphi_dR[4];
    double omega_adapt_raw, R_adapt_raw;
    double R_term_a, R_term_b;
    double theta_dot_w, theta_dot_R;
    double v_nu[4];
    double dot_psis_a, dot_psis_b;
    double dot_chi_a, dot_chi_b;
    double eaf_a, eaf_b;
    double psi_af_a, psi_af_b;
    double theta_hat;
    double phi_hat;
    double nu_pe;

    (void)s;

    if (!g_init) {
        observer_ipmsm_nlrs_init(NULL);
    }

    memset(out, 0, sizeof(*out));

    /* -----------------------------------------------------
     * 1) 输出 y = psi_sigma = Lq * i
     *    对应论文 y = Cx, C = [I 0]
     * ----------------------------------------------------- */
    y_psis_a = OBS_LQ * i_alpha;
    y_psis_b = OBS_LQ * i_beta;

    /* 输出误差 epsilon = psi_sigma - psi_sigma_hat */
    eps_a = y_psis_a - g_psis_hat[0];
    eps_b = y_psis_b - g_psis_hat[1];

    /* 低速时冻结符号，避免 sign(omega_hat)=0 导致位置公式抖动 */
    if (fabs(g_omega_hat_e) > OBS_SIGN_W_THR) {
        g_sign_omega = (g_omega_hat_e >= 0.0) ? 1.0 : -1.0;
    }

    /* -----------------------------------------------------
     * 2) 计算 ∂phi/∂theta(theta_hat)
     *    phi = [ -R i + omega J Lq i ; -omega J (u - R i) ]
     *
     *    ∂phi/∂omega = [ J Lq i ; -J(u - R i) ]
     *    ∂phi/∂R     = [ -i ; omega J i ]
     * ----------------------------------------------------- */
    J_mul(i_alpha, i_beta, &Ji_a, &Ji_b);

    /* u - R_hat i */
    {
        double tmp_a = u_alpha - g_R_hat * i_alpha;
        double tmp_b = u_beta  - g_R_hat * i_beta;
        J_mul(tmp_a, tmp_b, &JuRi_a, &JuRi_b);
    }

    dphi_domega[0] = OBS_LQ * Ji_a;
    dphi_domega[1] = OBS_LQ * Ji_b;
    dphi_domega[2] = -JuRi_a;
    dphi_domega[3] = -JuRi_b;

    dphi_dR[0] = -i_alpha;
    dphi_dR[1] = -i_beta;
    dphi_dR[2] = g_omega_hat_e * Ji_a;
    dphi_dR[3] = g_omega_hat_e * Ji_b;

    /* -----------------------------------------------------
     * 3) 参数更新律
     *    dot(theta_hat) = gamma * (nu_top^T * epsilon)
     *
     *    这里 nu_top 是 ν 的前 2 行，对应 C^T 的选择作用。
     * ----------------------------------------------------- */
    //omega_adapt_raw = -(g_nu[0][0] * eps_a + g_nu[1][0] * eps_b);
    //R_adapt_raw     = (g_nu[0][1] * eps_a + g_nu[1][1] * eps_b);

    omega_adapt_raw = -(g_nu[0][0] * eps_a + g_nu[1][0] * eps_b);

    R_term_a = g_nu[0][1] * eps_a;
    R_term_b = g_nu[1][1] * eps_b;
    R_adapt_raw = R_term_a + R_term_b;

    g_R_adapt_raw_f += Ts * OBS_R_ADAPT_LPF_W * (R_adapt_raw - g_R_adapt_raw_f);

    theta_dot_w = -OBS_GAMMA_W * omega_adapt_raw;
    theta_dot_R = OBS_GAMMA_R * R_adapt_raw;

    //禁用低通滤波电阻因子
    //theta_dot_R = OBS_GAMMA_R * g_R_adapt_raw_f;

    /* -----------------------------------------------------
     * 4) 状态观测器
     *    dot(x_hat) = Ac x_hat + phi(u,y,theta_hat)
     *               + lambda S^{-1} C^T epsilon + nu dot(theta_hat)
     *
     *    由于 S^{-1} C^T = [2I; I]^T，且 lambda=[lambda1 lambda2]^T，
     *    所以前两维反馈取 2*lambda1*epsilon，后两维取 lambda2*epsilon。
     * ----------------------------------------------------- */
    v_nu[0] = g_nu[0][0] * theta_dot_w + g_nu[0][1] * theta_dot_R;
    v_nu[1] = g_nu[1][0] * theta_dot_w + g_nu[1][1] * theta_dot_R;
    v_nu[2] = g_nu[2][0] * theta_dot_w + g_nu[2][1] * theta_dot_R;
    v_nu[3] = g_nu[3][0] * theta_dot_w + g_nu[3][1] * theta_dot_R;

    dot_psis_a = u_alpha - g_R_hat * i_alpha + g_chi_hat[0]
               + g_omega_hat_e * OBS_LQ * Ji_a
               + 2.0 * OBS_LAMBDA1 * eps_a
               + v_nu[0];
    dot_psis_b = u_beta  - g_R_hat * i_beta  + g_chi_hat[1]
               + g_omega_hat_e * OBS_LQ * Ji_b
               + 2.0 * OBS_LAMBDA1 * eps_b
               + v_nu[1];

    dot_chi_a = -g_omega_hat_e * JuRi_a
              + OBS_LAMBDA2 * eps_a
              + v_nu[2];
    dot_chi_b = -g_omega_hat_e * JuRi_b
              + OBS_LAMBDA2 * eps_b
              + v_nu[3];

    g_psis_hat[0] += Ts * dot_psis_a;
    g_psis_hat[1] += Ts * dot_psis_b;
    g_chi_hat[0]  += Ts * dot_chi_a;
    g_chi_hat[1]  += Ts * dot_chi_b;

    //g_omega_hat_e_obs += Ts * theta_dot_w;

    // if(g_debug_time < 0.1)
    //     g_R_hat       = OBS_RS_INIT;
    // else    
    //     g_R_hat       += Ts * theta_dot_R;
    //g_R_hat       = OBS_RS_INIT;
    //g_omega_hat_e = g_dbg_omega_true_e;

    g_R_hat       += Ts * theta_dot_R;
    g_omega_hat_e += Ts * theta_dot_w;

    g_omega_hat_e = clamp(g_omega_hat_e, -OBS_OMEGA_MAX, OBS_OMEGA_MAX);
    g_R_hat       = clamp(g_R_hat, OBS_RS_MIN, OBS_RS_MAX);

    /* -----------------------------------------------------
     * 5) 过滤回归器 ν 更新
     *    采用按块写开的形式：
     *      dot(nu_psi) = -2*lambda1*nu_psi + lambda1*nu_chi + dphi_psi/dtheta
     *      dot(nu_chi) = -lambda2*nu_psi + dphi_chi/dtheta
     * ----------------------------------------------------- */
    {
        int j;
        for (j = 0; j < 2; ++j) {
            double dphi0 = (j == 0) ? dphi_domega[0] : dphi_dR[0];
            double dphi1 = (j == 0) ? dphi_domega[1] : dphi_dR[1];
            double dphi2 = (j == 0) ? dphi_domega[2] : dphi_dR[2];
            double dphi3 = (j == 0) ? dphi_domega[3] : dphi_dR[3];

            double kcol = 1;//(j == 1) ? OBS_NU_R_SCALE : 1.0;
            double kup  = 1;//(j == 1) ? OBS_NU_R_UP_SCALE : 1.0;

            double dot_nu0 = -OBS_LAMBDA1 * g_nu[0][j] + g_nu[2][j] + dphi0;
            double dot_nu1 = -OBS_LAMBDA1 * g_nu[1][j] + g_nu[3][j] + dphi1;
            double dot_nu2 = -OBS_LAMBDA2 * g_nu[0][j] +              dphi2;
            double dot_nu3 = -OBS_LAMBDA2 * g_nu[1][j] +              dphi3;

            // double dot_nu0 = -2.0 * OBS_LAMBDA1 * g_nu[0][j] + OBS_LAMBDA1 * g_nu[2][j] + dphi0;
            // double dot_nu1 = -2.0 * OBS_LAMBDA1 * g_nu[1][j] + OBS_LAMBDA1 * g_nu[3][j] + dphi1;
            // double dot_nu2 = -OBS_LAMBDA2 * g_nu[0][j] + dphi2;
            // double dot_nu3 = -OBS_LAMBDA2 * g_nu[1][j] + dphi3;

            g_nu[0][j] += Ts * dot_nu0;
            g_nu[1][j] += Ts * dot_nu1;
            g_nu[2][j] += Ts * dot_nu2;
            g_nu[3][j] += Ts * dot_nu3;
        }
    }

    /* -----------------------------------------------------
     * 6) 用 chi_hat + omega_hat 恢复 active EMF，再转成 active-flux 方向
     *
     *    e_active_hat = -chi_hat - omega_hat * J * Lq i
     *    psi_active_dir ~ -sign(omega_hat) * J * e_active_hat
     *
     *    这样得到的方向与论文式(18)等价，但实现上更直观。
     * ----------------------------------------------------- */
    eaf_a = -g_chi_hat[0] - g_omega_hat_e * OBS_LQ * Ji_a;
    eaf_b = -g_chi_hat[1] - g_omega_hat_e * OBS_LQ * Ji_b;

    /* -sign(w) * J * e_active */
    psi_af_a =  g_sign_omega * eaf_b;
    psi_af_b = -g_sign_omega * eaf_a;

    phi_hat = sqrt(psi_af_a * psi_af_a + psi_af_b * psi_af_b);
    if (phi_hat < OBS_EPS) {
        phi_hat = OBS_EPS;
    }
    theta_hat = atan2(psi_af_b, psi_af_a);
    theta_hat = wrap_pm_pi(theta_hat);

    nu_pe = g_nu[0][0] * g_nu[0][0] + g_nu[1][0] * g_nu[1][0]
          + g_nu[0][1] * g_nu[0][1] + g_nu[1][1] * g_nu[1][1];

    /* -----------------------------------------------------
     * 7) 输出到现有工程 / logger
     * ----------------------------------------------------- */
    out->theta_e_hat = theta_hat;
    out->omega_e_hat = g_omega_hat_e;//g_omega_hat_e_obs;

    /* 这里把 active-flux 方向向量塞到 psi_hat，便于你直接画方向 */
    out->psi_alpha_hat = psi_af_a;
    out->psi_beta_hat  = psi_af_b;
    out->phi_hat       = phi_hat;

    out->R_hat = g_R_hat;
    out->beta_hat = g_omega_hat_e;
    out->detQ = nu_pe;

    /* 复用现有诊断槽位 */
    // out->eta1_hat = omega_adapt_raw;  /* 速度通道原始更新量，未乘 gamma */
    // out->eta2_hat = theta_dot_w;      /* 速度通道最终更新量，已乘 gamma */

    // out->q1 = eps_a;
    // out->q2 = eps_b;
    // out->q3 = g_psis_hat[0];
    // out->q4 = g_psis_hat[1];
    // out->q5 = g_chi_hat[0];
    // out->q6 = g_chi_hat[1];

    out->eta1_hat = R_term_a;         /* R通道 alpha 分量 */
    out->eta2_hat = R_term_b;         /* R通道 beta 分量 */

    out->q1 = eps_a;
    out->q2 = eps_b;
    
    out->q5 = g_nu[2][1];
    out->q6 = g_nu[3][1];

    out->yreg = phi_hat;

    out->z21 = g_nu[0][0];
    out->z22 = g_nu[1][0];
    out->xi1 = psi_af_a;
    out->xi2 = psi_af_b;

    /* eta1/eta2 对这篇论文没有意义，保留为 0 */
}
