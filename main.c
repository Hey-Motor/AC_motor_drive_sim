#include <math.h>
#include <stdio.h>

#include "config/params.h"
#include "config/config.h"
#include "control/foc_controller.h"
#include "control/im_observer_fo.h"
#include "io/adc_model.h"
#include "io/logger.h"
#include "plant/im_model.h"
#include "plant/inverter_avg.h"
#include "plant/pmsm_model.h"

double g_debug_time = 0.0;

typedef struct {
    int pass_speed;
    double speed_avg_err;
    double wmhat_rmse_last1s;
} CaseResult;

static CaseResult run_case(double speed_ref,
                           const char *csv_name,
                           int full_sensorless,
                           int observer_use_real_speed,
                           int print_progress)
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

    double id_ref = 0.0;
    double iq_ref = 0.0;
    double u_alpha_cmd = 0.0;
    double u_beta_cmd = 0.0;
    double u_alpha_act_next = 0.0;
    double u_beta_act_next = 0.0;
    double omega_m_avg_last_1s = 0.0;
    int omega_m_avg_count = 0;
    double wmhat_err_sq_sum = 0.0;

    CaseResult ret = {0, 0.0, 0.0};

    params_load_default(&p);
    im_init(&x);
    foc_init(&foc);

    im_observer_fo_init(&obs_state);
    im_observer_fo_set_speed_init(speed_ref * p.motor_im.pole_pairs);
    im_observer_fo_set_use_real_speed(observer_use_real_speed);

    u.u_alpha = 0.0;
    u.u_beta = 0.0;
    u.T_load = 0.0;

    im_get_output(&x, &u, &p.motor_im, &y);

    if (!logger_open(&lg, csv_name)) {
        printf("failed to open %s\n", csv_name);
        return ret;
    }

    {
        int total_steps = (int)(p.sim.t_end / p.sim.Ts_ctrl);
        int print_decim = (int)(0.2 / p.sim.Ts_ctrl);
        if (print_decim < 1) {
            print_decim = 1;
        }

        for (int k = 0; k < total_steps; ++k) {
            double t = k * p.sim.Ts_ctrl;
            double omega_m_hat;

            g_debug_time = t;
            u.T_load = 0.0;

            u.u_alpha = u_alpha_act_next;
            u.u_beta = u_beta_act_next;

            for (int n = 0; n < p.sim.substeps; ++n) {
                im_step_rk4(&x, &u, &p.motor_im, p.sim.dt_plant);
            }

            im_get_output(&x, &u, &p.motor_im, &y);
            adc_sample_currents(y.ia, y.ib, &p.adc, &adc);

            im_observer_fo_step(&obs_state, p.sim.Ts_ctrl, u.u_alpha, u.u_beta,
                                adc.i_alpha, adc.i_beta, &p.motor_im, y.omega_m, &obs);

            omega_m_hat = obs.omega_e_hat / p.motor_im.pole_pairs;

            if (print_progress && ((k % print_decim) == 0)) {
                printf("t=%.3f s, speed_ref=%.1f, speed_real=%.3f, speed_hat=%.3f\n",
                       t, speed_ref, y.omega_m, omega_m_hat);
            }

            if (full_sensorless) {
                foc_step(speed_ref, obs.theta_e_hat, omega_m_hat, &adc,
                         &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);
            } else {
                foc_step(speed_ref, obs.theta_e_hat, y.omega_m, &adc,
                         &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);
            }

            inverter_apply(u_alpha_cmd, u_beta_cmd,
                           adc.ia, adc.ib, -(adc.ia + adc.ib),
                           &p.inverter,
                           &u_alpha_act_next, &u_beta_act_next);

            if (t >= (p.sim.t_end - 1.0)) {
                double e = omega_m_hat - y.omega_m;
                omega_m_avg_last_1s += y.omega_m;
                wmhat_err_sq_sum += e * e;
                omega_m_avg_count += 1;
            }

            logger_write(&lg, t, speed_ref, id_ref, iq_ref,
                         u_alpha_cmd, u_beta_cmd, u_alpha_act_next, u_beta_act_next,
                         &y, &adc, &obs, 0.0, u.T_load);
        }
    }

    logger_close(&lg);

    if (omega_m_avg_count > 0) {
        double omega_avg = omega_m_avg_last_1s / (double)omega_m_avg_count;
        ret.speed_avg_err = omega_avg - speed_ref;
        ret.wmhat_rmse_last1s = sqrt(wmhat_err_sq_sum / (double)omega_m_avg_count);
        ret.pass_speed = (ret.speed_avg_err < 3.0 && ret.speed_avg_err > -3.0) ? 1 : 0;
    }

    return ret;
}

int main(void)
{
    const double speed_ref = 150.0;
    const char *csv_name = "run_im_single.csv";
    CaseResult r;

    im_observer_fo_set_adapt_gain_sched_enabled(0);
    foc_set_speed_gain_sched_enabled(0);
    im_observer_fo_set_adapt_err_lpf_enabled(0);

    printf("\n[Single Run] full sensorless, RK4 observer\n");
    printf("gain schedule: speed_loop=OFF, adapt_law=OFF, adapt_lpf=OFF\n");
    printf("save csv: %s\n", csv_name);

    r = run_case(speed_ref, csv_name, 1, 0, 1);
    printf("summary: speed_ref=%.1f, steady_speed_err=%.4f, wmhat_rmse=%.4f, result=%s\n",
           speed_ref, r.speed_avg_err, r.wmhat_rmse_last1s, r.pass_speed ? "PASS" : "FAIL");

    return 0;
}
