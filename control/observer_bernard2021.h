#ifndef OBSERVER_BERNARD2021_H
#define OBSERVER_BERNARD2021_H

#include "observer_select.h"

/* =========================================================
 * Bernard & Praly, TAC 2021
 * "Estimation of Position and Resistance of a Sensorless PMSM:
 *  A Nonlinear Luenberger Approach for a Nonobservable System"
 *
 * 这版代码按你当前工程接口做了一个“可集成的第一版”：
 * 1) 动态部分：严格按论文式(26)实现 3 组 lambda 滤波器；
 * 2) 静态部分：按论文 IV-A 的思路，在 1-D 电阻网格上搜索；
 * 3) 速度：额外实现论文式(38) 的速度估计器，给 FOC 提供 omega_hat；
 * 4) 为了不大改 core/types.h，内部状态仍采用 file-static 保存。
 *
 * 重要说明：
 * - 论文主体针对 nonsalient PMSM（SPMSM）。
 * - 因此本分支建议配合 USE_IPMSM=0 使用；
 * - 若在 IPMSM 上直接用，本实现只可视为“近似实验版”，
 *   不能宣称完全对应论文的严格模型。
 * ========================================================= */

void observer_bernard2021_init(ObserverState *s);
void observer_bernard2021_set_debug_truth(double theta_e_true, double omega_e_true);

void observer_bernard2021_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif
