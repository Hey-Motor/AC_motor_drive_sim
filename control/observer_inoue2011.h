#ifndef OBSERVER_INOUE2011_H
#define OBSERVER_INOUE2011_H

#include "observer_select.h"

/*
 * Inoue 2011 / Morimoto 2002 分支
 *
 * 接口保持和当前工程一致：
 *   - observer_inoue2011_init()
 *   - observer_inoue2011_set_debug_truth()
 *   - observer_inoue2011_step()
 *
 * 说明：
 * 1) 位置/速度无感本体：采用 Morimoto 2002 的 rotating gamma-delta
 *    extended-EMF least-order observer。
 * 2) 电阻在线辨识：采用 Inoue 2011 的 q-axis（在估计 gamma-delta 坐标下
 *    对应 delta 轴）RLS 标量辨识模型。
 * 3) 由于当前工程的 ObserverState 结构体是给其他 observer 复用的，
 *    本文件内部额外状态采用 file-static 方式保存，避免大改 core/types.h。
 */

void observer_inoue2011_init(ObserverState *s);
void observer_inoue2011_set_debug_truth(double theta_e_true, double omega_e_true);

void observer_inoue2011_step(
    ObserverState *s,
    double Ts,
    double u_alpha,
    double u_beta,
    double i_alpha,
    double i_beta,
    ObserverOutput *out
);

#endif
