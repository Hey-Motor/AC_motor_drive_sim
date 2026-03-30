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
                           int observer_use_real_speed)
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
    const double refs[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
    int first_fail = -1;

    im_observer_fo_set_adapt_gain_sched_enabled(1);
    foc_set_speed_gain_sched_enabled(1);
    im_observer_fo_set_adapt_err_lpf_enabled(0);

    printf("\n[Task1] Full sensorless + RK4 observer (10~90 rad/s)\n");
    for (int i = 0; i < (int)(sizeof(refs) / sizeof(refs[0])); ++i) {
        char name_buf[64];
        CaseResult r;
        snprintf(name_buf, sizeof(name_buf), "run_im_ref%.0f_full_sensorless_rk4.csv", refs[i]);
        r = run_case(refs[i], name_buf, 1, 0);
        printf("rk4_full_ref%.0f: speed_err=%.4f, rmse=%.4f, result=%s\n",
               refs[i], r.speed_avg_err, r.wmhat_rmse_last1s, r.pass_speed ? "PASS" : "FAIL");
        if ((!r.pass_speed) && first_fail < 0) {
            first_fail = i;
        }
    }
    if (first_fail < 0) {
        printf("rk4_full_scan: PASS up to 90 rad/s\n");
    } else {
        printf("rk4_full_scan: first FAIL at %.1f rad/s\n", refs[first_fail]);
    }

    {
        double rated_mech_ref = 2.0 * 3.14159265358979323846 * 50.0 / 2.0; /* rated electrical 50Hz, pole_pairs=2 */
        CaseResult r;
        printf("\n[Task2] With speed sensor, observer uses real speed (rated speed test)\n");
        r = run_case(rated_mech_ref,
                     "run_im_rated_speed_sensor_mode.csv",
                     0,  /* speed loop uses real speed */
                     1); /* observer uses real speed */
        printf("rated_sensor_mode: ref=%.4f, speed_err=%.4f, rmse=%.4f, result=%s\n",
               rated_mech_ref, r.speed_avg_err, r.wmhat_rmse_last1s, r.pass_speed ? "PASS" : "FAIL");
    }

    {
        typedef struct {
            const char *tag;
            int adapt_sched_on;
            int speed_sched_on;
            int adapt_lpf_on;
        } AblationCfg;

        const AblationCfg cfgs[] = {
            {"none",               0, 0, 0},
            {"adapt_sched_only",   1, 0, 0},
            {"speed_sched_only",   0, 1, 0},
            {"adapt_lpf_only",     0, 0, 1},
            {"all_three",          1, 1, 1},
            {"current_default",    1, 1, 0}
        };

        printf("\n[Task3] 150 rad/s full-sensorless ablation\n");
        printf("pass criterion: steady-state speed error in last 1s within +/-3 rad/s\n");

        for (int i = 0; i < (int)(sizeof(cfgs) / sizeof(cfgs[0])); ++i) {
            char name_buf[96];
            CaseResult r;

            im_observer_fo_set_adapt_gain_sched_enabled(cfgs[i].adapt_sched_on);
            foc_set_speed_gain_sched_enabled(cfgs[i].speed_sched_on);
            im_observer_fo_set_adapt_err_lpf_enabled(cfgs[i].adapt_lpf_on);

            snprintf(name_buf, sizeof(name_buf), "run_im_150_ablation_%s.csv", cfgs[i].tag);
            r = run_case(150.0, name_buf, 1, 0);

            printf("ablation %-16s | adapt_sched=%d speed_sched=%d adapt_lpf=%d"
                   " | speed_err=%.4f rmse=%.4f | pass150=%s\n",
                   cfgs[i].tag,
                   cfgs[i].adapt_sched_on,
                   cfgs[i].speed_sched_on,
                   cfgs[i].adapt_lpf_on,
                   r.speed_avg_err,
                   r.wmhat_rmse_last1s,
                   r.pass_speed ? "YES" : "NO");
        }
    }

    {
        const double scan_start = 100.0;
        const double scan_end = 300.0;
        const double scan_step = 10.0;
        double max_pass_ref = 0.0;
        int has_fail = 0;

        printf("\n[Task4] Full-sensorless max-speed scan with current default setup\n");
        im_observer_fo_set_adapt_gain_sched_enabled(1);
        foc_set_speed_gain_sched_enabled(1);
        im_observer_fo_set_adapt_err_lpf_enabled(0);

        for (double ref = scan_start; ref <= scan_end + 1e-9; ref += scan_step) {
            char name_buf[96];
            CaseResult r;
            snprintf(name_buf, sizeof(name_buf), "run_im_maxscan_ref%.0f.csv", ref);
            r = run_case(ref, name_buf, 1, 0);

            printf("maxscan ref=%.0f | speed_err=%.4f rmse=%.4f | result=%s\n",
                   ref, r.speed_avg_err, r.wmhat_rmse_last1s, r.pass_speed ? "PASS" : "FAIL");

            if (r.pass_speed) {
                max_pass_ref = ref;
            } else {
                has_fail = 1;
                break;
            }
        }

        if (has_fail) {
            printf("maxscan summary: max PASS reference = %.0f rad/s\n", max_pass_ref);
        } else {
            printf("maxscan summary: PASS up to %.0f rad/s (scan limit)\n", scan_end);
        }
    }

    return 0;
}
