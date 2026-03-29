#include <math.h>
#include <stdio.h>

#include "config/params.h"
#include "control/foc_controller.h"
#include "control/observer_select.h"
#include "io/adc_model.h"
#include "io/logger.h"
#include "plant/inverter_avg.h"
#include "plant/pmsm_model.h"

/* 仅供 FOC 调试注入使用 */
double g_debug_time = 0.0;


/* 角度误差统一包到 [-pi, pi] */
static double wrap_pm_pi(double x)
{
    while (x > 3.14159265358979323846) {
        x -= 2.0 * 3.14159265358979323846;
    }
    while (x < -3.14159265358979323846) {
        x += 2.0 * 3.14159265358979323846;
    }
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

    double speed_ref = 50.0;
    double id_ref = 0.0;
    double iq_ref = 0.0;

    /* 控制器算出的电压命令 */
    double u_alpha_cmd = 0.0;
    double u_beta_cmd = 0.0;

    /* 下一拍真正施加到电机的平均电压 */
    double u_alpha_act_next = 0.0;
    double u_beta_act_next = 0.0;

    /* 收敛评估统计（t >= 1.0s） */
    double max_abs_theta_err_e_post1s = 0.0;
    double max_abs_omega_err_m_post1s = 0.0;

    params_load_default(&p);
    pmsm_init(&x);
    foc_init(&foc);
    observer_init(&obs_state);

    u.u_alpha = 0.0;
    u.u_beta = 0.0;
    u.T_load = 0.0;

    pmsm_get_output(&x, &u, &p.motor, &y);

    if (!logger_open(&lg, "run.csv")) {
        printf("failed to open run.csv\n");
        return 1;
    }

    printf("t,omega_m,omega_e,omega_m_hat,omega_e_hat,theta_err_e,omega_err_m\n");

    {
        int total_steps = (int)(p.sim.t_end / p.sim.Ts_ctrl);

        for (int k = 0; k < total_steps; ++k) {
            double t = k * p.sim.Ts_ctrl;
            double omega_m_hat;
            double theta_err_e;
            double omega_err_m;

            g_debug_time = t;

            if(t>3)
                u.T_load = 2.0;

            /* 固定速度参考：50 rad/s (mechanical) */
            speed_ref = 20.0;

            /*
             * 按要求注释掉 Rs 随时间变化逻辑：
             * if (t < 5.0) { p.motor.Rs = 1.025; }
             * else if (t < 10.0) { p.motor.Rs = 2.5; }
             * else { p.motor.Rs = 3.5; }
             */

            /* 单采单更时序：本拍开始施加上一拍锁存的实际电压 */
            u.u_alpha = u_alpha_act_next;
            u.u_beta = u_beta_act_next;

            for (int n = 0; n < p.sim.substeps; ++n) {
                pmsm_step_rk4(&x, &u, &p.motor, p.sim.dt_plant);
            }

            pmsm_get_output(&x, &u, &p.motor, &y);

            adc_sample_currents(y.ia, y.ib, &p.adc, &adc);

            /*
             * 本阶段说明：
             * 1) 只实现 Piippo2008 基础 adaptive observer；
             * 2) HF injection 后续再做；
             * 3) 当前 FOC 由传感器/真值驱动，observer 仅用于独立评估。
             */
            observer_set_debug_truth(y.theta_e, y.omega_e);
            observer_step(&obs_state, p.sim.Ts_ctrl, u.u_alpha, u.u_beta,
                          adc.i_alpha, adc.i_beta, &obs);

            omega_m_hat = obs.omega_e_hat / p.motor.pole_pairs;

            /* 0.2s 前用真实速度/位置，0.2s 后切换到 observer 闭环 */
            if (t < 0.2) {
                foc_step(speed_ref, y.theta_e, y.omega_m, &adc,
                         &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);
            } else {
                foc_step(speed_ref, obs.theta_e_hat, omega_m_hat, &adc,
                         &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);
            }
            
            inverter_apply(u_alpha_cmd, u_beta_cmd,
                           adc.ia, adc.ib, -(adc.ia + adc.ib),
                           &p.inverter,
                           &u_alpha_act_next, &u_beta_act_next);

            theta_err_e = wrap_pm_pi(obs.theta_e_hat - y.theta_e);
            omega_err_m = omega_m_hat - y.omega_m;

            if (t >= 1.0) {
                if (fabs(theta_err_e) > max_abs_theta_err_e_post1s) {
                    max_abs_theta_err_e_post1s = fabs(theta_err_e);
                }
                if (fabs(omega_err_m) > max_abs_omega_err_m_post1s) {
                    max_abs_omega_err_m_post1s = fabs(omega_err_m);
                }
            }

            logger_write(&lg, t, speed_ref, id_ref, iq_ref,
                         u_alpha_cmd, u_beta_cmd, u_alpha_act_next, u_beta_act_next,
                         &y, &adc, &obs, theta_err_e, u.T_load);

            if ((k % 100) == 0) {
                printf("t=%.4f,wm=%.4f,we=%.4f,wm_hat=%.4f,we_hat=%.4f,theta_err_e=%.4f,omega_err_m=%.4f\n",
                       t, y.omega_m, y.omega_e, omega_m_hat, obs.omega_e_hat,
                       theta_err_e, omega_err_m);
            }
        }
    }

    logger_close(&lg);

    printf("post_1s_max_abs_omega_err_m=%.6f (threshold 5.0)\n", max_abs_omega_err_m_post1s);
    printf("post_1s_max_abs_theta_err_e=%.6f (threshold 0.5)\n", max_abs_theta_err_e_post1s);

    if ((max_abs_omega_err_m_post1s < 5.0) && (max_abs_theta_err_e_post1s < 0.5)) {
        printf("convergence_check=PASS\n");
    } else {
        printf("convergence_check=FAIL\n");
    }

    return 0;
}
