/* =========================================================
 * observer_select.c
 *
 * 作用：
 * 统一观测器公共接口。
 * main.c / FOC / logger 永远只调用：
 *   observer_init()
 *   observer_set_debug_truth()
 *   observer_step()
 *
 * 真正使用哪个观测器，由 config.h 中的
 *   OBSERVER_MODE_DEFAULT
 * 决定。
 *
 * 以后新增观测器，只需要：
 * 1) 新建 observer_xxx.c / observer_xxx.h
 * 2) 在 config.h 里增加宏
 * 3) 在本文件里加一个分支
 * ========================================================= */

#include "../config/config.h"
#include "observer_select.h"
#include "observer_fluxphi.h"
#include "observer_pebo_drem.h"
#include "observer_bobtsov2015.h"
#include "observer_ipmsm_nlrs.h"
#include "observer_inoue2011.h"
#include "observer_bernard2021.h"
#include "observer_eta_only.h"
#include "observer_piippo2008.h"

/* =========================================================
 * 观测器公共接口
 * 你主循环 main.c / FOC / logger 都继续只认这 3 个函数
 * 具体用哪个观测器，由宏 OBSERVER_MODE_DEFAULT 决定
 * ========================================================= */

void observer_init(ObserverState *s)
{
#if (OBSERVER_MODE_DEFAULT == OBSERVER_FLUXPHI)
    observer_fluxphi_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PEBO_DREM)
    observer_pebo_drem_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BOBTSOV2015)
    observer_bobtsov2015_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_IPMSM_NLRS)
    observer_ipmsm_nlrs_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_INOUE2011)
    observer_inoue2011_init(s);        
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BERNARD2021)
    observer_bernard2021_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_ETA_ONLY)
    observer_eta_only_init(s);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PIIPPO2008)
    observer_piippo2008_init(s);
#else
#error "Unknown observer mode"
#endif
}

void observer_set_debug_truth(double theta_e_true, double omega_e_true)
{
#if (OBSERVER_MODE_DEFAULT == OBSERVER_FLUXPHI)
    observer_fluxphi_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PEBO_DREM)
    observer_pebo_drem_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BOBTSOV2015)
    observer_bobtsov2015_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_IPMSM_NLRS)
    observer_ipmsm_nlrs_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_INOUE2011)
    observer_inoue2011_set_debug_truth(theta_e_true, omega_e_true);        
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BERNARD2021)
    observer_bernard2021_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_ETA_ONLY)
    observer_eta_only_set_debug_truth(theta_e_true, omega_e_true);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PIIPPO2008)
    observer_piippo2008_set_debug_truth(theta_e_true, omega_e_true);
#else
#error "Unknown observer mode"
#endif
}

void observer_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
)
{
#if (OBSERVER_MODE_DEFAULT == OBSERVER_FLUXPHI)
    observer_fluxphi_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PEBO_DREM)
    observer_pebo_drem_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BOBTSOV2015)
    observer_bobtsov2015_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_IPMSM_NLRS)
    observer_ipmsm_nlrs_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_INOUE2011)
    observer_inoue2011_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);            
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_BERNARD2021)
    observer_bernard2021_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_ETA_ONLY)
    observer_eta_only_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#elif (OBSERVER_MODE_DEFAULT == OBSERVER_PIIPPO2008)
    observer_piippo2008_step(s, Ts, u_alpha, u_beta, i_alpha, i_beta, out);
#else
#error "Unknown observer mode"
#endif
}
