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
    double kp_r;
    double ki_r;
} RsTuneConfig;

typedef struct {
    int pass_speed;
    int pass_rs;
    double speed_avg_err;
    double wmhat_rmse_last1s;
    double rs_avg_last1s;
    double rs_err_last1s;
    double rs_track_mae;
    double t_conv_after_2s;
    double t_conv_after_4s;
} CaseResult;

static double rs_profile(double t, double rs_nom)
{
    if (t < 2.0) return rs_nom;
    if (t < 6.0) return 1.5 * rs_nom;
    return 0.5 * rs_nom;
}

static CaseResult run_case(double speed_ref,
                           const RsTuneConfig *cfg,
                           const char *csv_name,
                           double t_end,
                           int save_csv,
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
    double rs_hat_avg_last_1s = 0.0;
    double rs_abs_err_sum = 0.0;
    int rs_err_count = 0;
    int avg_count = 0;
    double wmhat_err_sq_sum = 0.0;
    double rs_nom = 0.0;
    double conv_hold_2s = 0.0;
    double conv_hold_4s = 0.0;

    CaseResult ret = {0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0};

    params_load_default(&p);
    rs_nom = p.motor_im.Rs;
    p.sim.t_end = t_end;
    im_init(&x);
    foc_init(&foc);

    im_observer_fo_init(&obs_state);
    im_observer_fo_set_speed_init(speed_ref * p.motor_im.pole_pairs);
    im_observer_fo_set_use_real_speed(0);
    im_observer_fo_set_adapt_gain_sched_enabled(0);
    im_observer_fo_set_adapt_err_lpf_enabled(0);
    im_observer_fo_set_rs_adapt_enabled(1);
    im_observer_fo_set_rs_adapt_gains(cfg->kp_r, cfg->ki_r);

    foc_set_speed_gain_sched_enabled(0);
    u.u_alpha = 0.0;
    u.u_beta = 0.0;
    u.T_load = 0.0;

    im_get_output(&x, &u, &p.motor_im, &y);

    if (save_csv) {
        if (!logger_open(&lg, csv_name)) {
            printf("failed to open %s\n", csv_name);
            return ret;
        }
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
            double rs_true;

            g_debug_time = t;
            u.T_load = 0.0;
            rs_true = rs_profile(t, rs_nom);
            /* 若不想启用变Rs工况，可改为：rs_true = rs_nom; */
            p.motor_im.Rs = rs_true;

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
                printf("t=%.3f s, speed_ref=%.1f, speed_real=%.3f, speed_hat=%.3f, Rs_true=%.4f, Rs_hat=%.4f\n",
                       t, speed_ref, y.omega_m, omega_m_hat, rs_true, obs.R_hat);
            }

            foc_step(speed_ref, obs.theta_e_hat, omega_m_hat, &adc,
                     &p, &foc, &u_alpha_cmd, &u_beta_cmd, &id_ref, &iq_ref);

            inverter_apply(u_alpha_cmd, u_beta_cmd,
                           adc.ia, adc.ib, -(adc.ia + adc.ib),
                           &p.inverter,
                           &u_alpha_act_next, &u_beta_act_next);

            if (t >= (p.sim.t_end - 1.0)) {
                double e = omega_m_hat - y.omega_m;
                omega_m_avg_last_1s += y.omega_m;
                rs_hat_avg_last_1s += obs.R_hat;
                wmhat_err_sq_sum += e * e;
                avg_count += 1;
            }
            if (t >= 2.0) {
                double rs_abs_err = fabs(obs.R_hat - rs_true);
                rs_abs_err_sum += rs_abs_err;
                rs_err_count += 1;

                if (t >= 2.0 && t < 6.0 && ret.t_conv_after_2s < 0.0) {
                    conv_hold_2s = (rs_abs_err <= 0.08) ? (conv_hold_2s + p.sim.Ts_ctrl) : 0.0;
                    if (conv_hold_2s >= 0.2) {
                        ret.t_conv_after_2s = t - conv_hold_2s + p.sim.Ts_ctrl;
                    }
                } else if (t >= 6.0 && ret.t_conv_after_4s < 0.0) {
                    conv_hold_4s = (rs_abs_err <= 0.08) ? (conv_hold_4s + p.sim.Ts_ctrl) : 0.0;
                    if (conv_hold_4s >= 0.2) {
                        ret.t_conv_after_4s = t - conv_hold_4s + p.sim.Ts_ctrl;
                    }
                }
            } 

            if (save_csv) {
                logger_write(&lg, t, speed_ref, id_ref, iq_ref,
                             u_alpha_cmd, u_beta_cmd, u_alpha_act_next, u_beta_act_next,
                             &y, &adc, &obs, 0.0, u.T_load);
            }
        }
    }

    if (save_csv) {
        logger_close(&lg);
    }

    if (avg_count > 0) {
        double omega_avg = omega_m_avg_last_1s / (double)avg_count;
        double rs_avg = rs_hat_avg_last_1s / (double)avg_count;
        ret.speed_avg_err = omega_avg - speed_ref;
        ret.wmhat_rmse_last1s = sqrt(wmhat_err_sq_sum / (double)avg_count);
        ret.rs_avg_last1s = rs_avg;
        ret.rs_err_last1s = rs_avg - p.motor_im.Rs;
        ret.rs_track_mae = (rs_err_count > 0) ? (rs_abs_err_sum / (double)rs_err_count) : 0.0;

        ret.pass_speed = (ret.speed_avg_err < 3.0 && ret.speed_avg_err > -3.0) ? 1 : 0;
        ret.pass_rs = (ret.rs_track_mae <= 0.08) ? 1 : 0;
    }

    return ret;
}

int main(void)
{
    const double speed_ref = 150.0;
    const char *csv_name = "run_im_var_rs.csv";
    const double kp_list[] = {0.02, 0.05, 0.10, 0.20, 0.50};
    const double ki_list[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    const int nkp = (int)(sizeof(kp_list) / sizeof(kp_list[0]));
    const int nki = (int)(sizeof(ki_list) / sizeof(ki_list[0]));
    RsTuneConfig best_cfg = {kp_list[0], ki_list[0]};
    CaseResult best = {0, 0, 0.0, 0.0, 0.0, 1e9};

    printf("\n[Var Rs Sweep] full sensorless DFOC (RK4 observer with Rs_hat in observer model)\n");
    printf("criterion: Rs tracking MAE (t>=2s) <= 0.08 ohm\n");

    for (int i = 0; i < nkp; ++i) {
        for (int j = 0; j < nki; ++j) {
            RsTuneConfig cfg = {kp_list[i], ki_list[j]};
            CaseResult r = run_case(speed_ref, &cfg, 0, 2.0, 0, 0);
            printf("sweep: kp_r=%.3f, ki_r=%.3f -> Rs_hat=%.4f, Rs_err=%.4f, speed_err=%.4f\n",
                   cfg.kp_r, cfg.ki_r, r.rs_avg_last1s, r.rs_err_last1s, r.speed_avg_err);
            if (fabs(r.rs_err_last1s) < fabs(best.rs_err_last1s)) {
                best = r;
                best_cfg = cfg;
            }
        }
    }

    printf("\n[Best Rs adaptation config] kp_r=%.3f, ki_r=%.3f\n",
           best_cfg.kp_r, best_cfg.ki_r);

    {
        CaseResult final_result = run_case(speed_ref, &best_cfg, csv_name, 10.0, 1, 1);
        printf("var-Rs summary: speed_ref=%.1f, speed_err=%.4f, wmhat_rmse=%.4f, Rs_hat_last1s=%.4f, Rs_err_last1s=%.4f, Rs_track_mae=%.4f, speed=%s, Rs=%s\n",
               speed_ref,
               final_result.speed_avg_err,
               final_result.wmhat_rmse_last1s,
               final_result.rs_avg_last1s,
               final_result.rs_err_last1s,
               final_result.rs_track_mae,
               final_result.pass_speed ? "PASS" : "FAIL",
               final_result.pass_rs ? "PASS" : "FAIL");
        if (final_result.t_conv_after_2s >= 0.0) {
            printf("Rs convergence after 2s switch: converged at t=%.3f s\n", final_result.t_conv_after_2s);
        } else {
            printf("Rs convergence after 2s switch: NOT converged within [2,6)s\n");
        }
        if (final_result.t_conv_after_4s >= 0.0) {
            printf("Rs convergence after 6s switch: converged at t=%.3f s\n", final_result.t_conv_after_4s);
        } else {
            printf("Rs convergence after 6s switch: NOT converged within [6,10]s\n");
        }
        printf("save csv: %s\n", csv_name);
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
