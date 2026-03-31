#include <stdio.h>
#include <math.h>

#include "config/params.h"
#include "plant/im_model.h"
#include "plant/inverter_avg.h"
#include "control/foc_controller.h"
#include "control/im_observer_fo.h"
#include "io/adc_model.h"
#include "io/logger.h"

/* =========================
 * 单文件单入口版本（仅 1 个 main，无 run_case）
 * ========================= */
/* 仅供 FOC 调试注入使用 */
double g_debug_time = 0.0;

/* 角度误差统一包到 [-pi, pi] */
static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) x -= 2.0 * 3.14159265358979323846;
    while (x < -3.14159265358979323846) x += 2.0 * 3.14159265358979323846;
    return x;
}

int main(void)
{
    Params p;
    PlantState x;
    PlantInput u;
    PlantOutput y;

    FocState foc;
    ObserverState obs_state;
    ObserverOutput obs;
    AdcSample adc;
    CsvLogger lg;

    /* 工况：IM 无速度传感器闭环（速度/角度都用观测器） */
    double speed_ref = 150.0;
    double id_ref = 0.0, iq_ref = 0.0;

    /* 控制器算出的电压命令 */
    double u_alpha_cmd = 0.0, u_beta_cmd = 0.0;

    /* 下一拍真正施加到电机的平均电压 */
    double u_alpha_act_next = 0.0, u_beta_act_next = 0.0;

    params_load_default(&p);
    p.sim.t_end = 10.0;

    im_init(&x);
    foc_init(&foc);
    im_observer_fo_init(&obs_state);

    /* 观测器设置 */
    im_observer_fo_set_speed_init(speed_ref * p.motor_im.pole_pairs);
    im_observer_fo_set_use_real_speed(0);          /* 观测器内部使用估计速度 */
    im_observer_fo_set_adapt_gain_sched_enabled(0);
    im_observer_fo_set_adapt_err_lpf_enabled(0);
    im_observer_fo_set_rs_adapt_enabled(1);
    im_observer_fo_set_rs_adapt_gains(0.02, 0.5);

    u.u_alpha = 0.0;
    u.u_beta  = 0.0;
    u.T_load  = 0.0;

    im_get_output(&x, &u, &p.motor_im, &y);

    if (!logger_open(&lg, "run_im_simple.csv")) {
        printf("failed to open run_im_simple.csv\n");
        return 1;
    }

    {
        int total_steps = (int)(p.sim.t_end / p.sim.Ts_ctrl);

        for (int k = 0; k < total_steps; ++k) {
            double t = k * p.sim.Ts_ctrl;
            double omega_m_hat;
            double theta_err;

            g_debug_time = t;

            /* 变 Rs 工况（若不需要，直接把 rs_true 固定为 p.motor_im.Rs 即可） */
            if (t < 2.0) {
                p.motor_im.Rs = 0.567;
            } else if (t < 6.0) {
                p.motor_im.Rs = 1.5 * 0.567;
            } else {
                p.motor_im.Rs = 0.5 * 0.567;
            }

            /* 本拍开始时，施加上一拍实际平均电压 */
            u.u_alpha = u_alpha_act_next;
            u.u_beta  = u_beta_act_next;

            /* 一个控制周期内连续对象积分 */
            for (int n = 0; n < p.sim.substeps; ++n) {
                im_step_rk4(&x, &u, &p.motor_im, p.sim.dt_plant);
            }

            /* 周期末采样 */
            im_get_output(&x, &u, &p.motor_im, &y);
            adc_sample_currents(y.ia, y.ib, &p.adc, &adc);

            /* IM 全阶磁链观测器 */
            im_observer_fo_step(&obs_state, p.sim.Ts_ctrl,
                                u.u_alpha, u.u_beta,
                                adc.i_alpha, adc.i_beta,
                                &p.motor_im, y.omega_m,
                                &obs);

            omega_m_hat = obs.omega_e_hat / p.motor_im.pole_pairs;

            /* 无速度传感器闭环 */
            foc_step(speed_ref, obs.theta_e_hat, omega_m_hat, &adc,
                     &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);

            /* 逆变器平均模型：得到下一拍实际电压 */
            inverter_apply(u_alpha_cmd, u_beta_cmd,
                           adc.ia, adc.ib, -(adc.ia + adc.ib),
                           &p.inverter,
                           &u_alpha_act_next, &u_beta_act_next);

            theta_err = wrap_pm_pi(y.theta_e - obs.theta_e_hat);

            logger_write(&lg, t, speed_ref, id_ref, iq_ref,
                         u_alpha_cmd, u_beta_cmd, u_alpha_act_next, u_beta_act_next,
                         &y, &adc, &obs, theta_err, u.T_load);

            if ((k % 100) == 0) {
                printf("t=%.4f, Rs_real=%.4f, Rs_hat=%.4f, wm=%.4f, wm_hat=%.4f, theta_err=%.6f\n",
                       t, p.motor_im.Rs, obs.R_hat, y.omega_m, omega_m_hat, theta_err);
            }
        }
    }

    logger_close(&lg);
    return 0;
}
