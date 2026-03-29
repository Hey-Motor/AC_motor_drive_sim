#ifndef OBSERVER_SELECT_H
#define OBSERVER_SELECT_H

#include "../core/types.h"

/*
 * 观测器接口说明
 *
 * 你将来改真实观测器时，只需要：
 * 1) 保持这个头文件的函数接口不变
 * 2) 修改 observer_dummy.c 里的实现
 *
 * 这样 main.c / foc_controller.c / logger.c 都不用改
 */

void observer_init(ObserverState *s);

/* 当前仅供模板调试使用：
 * 现在 dummy 观测器用真值透传，所以需要从 main.c 灌入真值
 * 以后换成真实观测器后，这个函数可以删掉或留空
 */
void observer_set_debug_truth(double theta_e_true, double omega_e_true);

/*
 * 输入:
 *   u_alpha, u_beta : 本拍施加给电机的 alpha-beta 电压
 *   i_alpha, i_beta : 本拍采样到的 alpha-beta 电流
 *
 * 输出:
 *   theta_e_hat, omega_e_hat
 */
void observer_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif