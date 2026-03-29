#ifndef INVERTER_AVG_H
#define INVERTER_AVG_H

#include "../core/types.h"

/*
 * 平均逆变器模型
 *
 * 输入:
 *   u_alpha_cmd, u_beta_cmd : 控制器给出的电压指令
 *   ia, ib, ic              : 当前采样到的三相电流（用于非线性平均误差）
 *   inv                     : 逆变器参数
 *
 * 输出:
 *   u_alpha_act, u_beta_act : 下一拍实际施加到电机的平均电压
 */
void inverter_apply(
    double u_alpha_cmd,
    double u_beta_cmd,
    double ia,
    double ib,
    double ic,
    const InverterParams *inv,
    double *u_alpha_act,
    double *u_beta_act
);

#endif