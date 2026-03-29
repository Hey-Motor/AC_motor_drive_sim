#include "observer_eta_only.h"
#include <string.h>

/*
 * eta-only observer skeleton
 * 说明：当前文件仅用于后续实验接线。
 */

/* 仅供接口占位：记录调试真值，当前骨架不使用真实算法。 */
static double g_dbg_theta_e_true = 0.0;
static double g_dbg_omega_e_true = 0.0;

void observer_eta_only_init(ObserverState *s)
{
    if (s != NULL) {
        memset(s, 0, sizeof(*s));
        s->initialized = 1;
    }
}

void observer_eta_only_set_debug_truth(double theta_e_true, double omega_e_true)
{
    g_dbg_theta_e_true = theta_e_true;
    g_dbg_omega_e_true = omega_e_true;
}

void observer_eta_only_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
    (void)Ts;
    (void)u_alpha;
    (void)u_beta;
    (void)i_alpha;
    (void)i_beta;

    if ((s == NULL) || (out == NULL)) {
        return;
    }
    if (!s->initialized) {
        observer_eta_only_init(s);
    }

    memset(out, 0, sizeof(*out));

    /*
     * 注意：这里只是骨架，尚未实现真实算法。
     * 当前仅输出占位值，便于先完成 observer_select 接线与编译联通。
     */
    out->theta_e_hat = g_dbg_theta_e_true;
    out->omega_e_hat = g_dbg_omega_e_true;
}
