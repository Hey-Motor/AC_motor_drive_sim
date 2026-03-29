#include "inverter_avg.h"
#include "../control/transforms.h"

/* 把数值限制在 [-1, 1]，用于构造简化非线性曲线 */
static double clamp1(double x)
{
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

void inverter_apply(
    double u_alpha_cmd,
    double u_beta_cmd,
    double ia,
    double ib,
    double ic,
    const InverterParams *inv,
    double *u_alpha_act,
    double *u_beta_act
)
{
    double alpha = u_alpha_cmd;
    double beta  = u_beta_cmd;

    /* 两电平逆变器在线性区可输出的 alpha-beta 最大电压幅值 */
    double vmax = inv->Vdc / 1.7320508075688772;

    /* ============= 模式 0：理想直通 =============
     * 不考虑母线限制，不考虑非线性
     * 主要用于最理想的算法基线验证
     */
    if (inv->mode == 0) {
        *u_alpha_act = alpha;
        *u_beta_act  = beta;
        return;
    }

    /* ============= 模式 1 / 2：先做平均模型圆限幅 =============
     * 这是最基本的“平均逆变器”
     */
    circle_limit(&alpha, &beta, vmax);

    /* ============= 模式 1：平均模型 + 圆限幅 =============
     * 不加入电流方向相关非线性
     */
    if (inv->mode == 1 || !inv->enable_nl) {
        *u_alpha_act = alpha;
        *u_beta_act  = beta;
        return;
    }

    /* ============= 模式 2：平均模型 + 非线性平均误差 =============
     * 做法：
     * 1) 先把 alpha-beta 电压变成三相相电压
     * 2) 对每相引入与相电流方向相关的平均误差
     * 3) 去掉零序（共模）分量
     * 4) 再转回 alpha-beta
     *
     * 这里的“去零序”是因为三相同时加同一个偏置，只会形成共模，
     * 对当前不接中性点的三相电机主要响应的线电压差影响不大。
     */
    {
        double ua, ub, uc;
        double verr_sat = inv->nl_sat_ratio * inv->Vdc;
        double ith = (inv->nl_current_thr > 1e-9) ? inv->nl_current_thr : 1.0;
        double uavg;

        inv_clarke(alpha, beta, &ua, &ub, &uc);

        /* 简化版梯形/饱和平均误差
         * 小电流区近似线性，大电流区饱和
         */
        ua -= verr_sat * clamp1(ia / ith);
        ub -= verr_sat * clamp1(ib / ith);
        uc -= verr_sat * clamp1(ic / ith);

        /* 去零序 */
        uavg = (ua + ub + uc) / 3.0;
        ua -= uavg;
        ub -= uavg;
        uc -= uavg;

        /* 转回 alpha-beta */
        alpha = ua;
        beta  = (ub - uc) / 1.7320508075688772;

        /* 再限一次幅值，避免非线性误差后超母线能力 */
        circle_limit(&alpha, &beta, vmax);
    }

    *u_alpha_act = alpha;
    *u_beta_act  = beta;
}