#include "params.h"
#include "config.h"

void params_load_default(Params *p)
{
    p->sim.Ts_ctrl   = 1.0 / CTRL_FREQ_HZ;
    p->sim.substeps  = PLANT_SUBSTEPS;
    p->sim.dt_plant  = p->sim.Ts_ctrl / p->sim.substeps;
    p->sim.t_end     = 5.0;
    p->sim.log_decim = LOG_DECIM;

    p->motor.pole_pairs = 3;
    p->motor.Rs    = 1.025;
#if USE_IPMSM
    p->motor.Ld    = 9e-3;
    p->motor.Lq    = 18.5e-3;
#else
    p->motor.Ld    = 9e-3;
    p->motor.Lq    = 9e-3;
#endif
    p->motor.psi_f = 0.249;
    p->motor.J     = 0.002379;
    p->motor.B     = 1e-4;

    /* IM 参数（图示额定参数） */
    p->motor_im.pole_pairs = 2;
    p->motor_im.Rs = 0.567;
    p->motor_im.Rr = 0.441;
    p->motor_im.Lm = 110.1e-3;
    p->motor_im.Ls = 114.1e-3;
    p->motor_im.Lr = 114.1e-3;
    p->motor_im.J  = 0.02;
    p->motor_im.B  = 1e-4;

    p->inverter.Vdc            = 540.0;
    p->inverter.nl_sat_ratio   = 0.05;   /* 默认 2% Vdc */
    p->inverter.nl_current_thr = 2.0;    /* 电流阈值 2A */
    p->inverter.mode           = INVERTER_MODE_DEFAULT;
    p->inverter.enable_nl      = USE_INVERTER_NL;

    p->adc.ia_offset = 0.0;
    p->adc.ib_offset = 0.0;
    p->adc.adc_lsb   = 0.0;
    p->adc.enable_offset = 0;
    p->adc.enable_quant  = 0;

    p->ctrl.kp_id = 8.0;
    p->ctrl.ki_id = 800.0;
    p->ctrl.kp_iq = 8.0;
    p->ctrl.ki_iq = 800.0;
    p->ctrl.kp_spd = 0.1;
    p->ctrl.ki_spd = 2.0;
    p->ctrl.iq_limit = 10.0;
    p->ctrl.vdq_limit = 520.0;

    /* IM DFOC 控制器参数（首版调试值） */
    p->ctrl_im.kp_id = 8.0;
    p->ctrl_im.ki_id = 400.0;
    p->ctrl_im.kp_iq = 8.0;
    p->ctrl_im.ki_iq = 400.0;
    p->ctrl_im.kp_spd = 0.2;
    p->ctrl_im.ki_spd = 4.0;
    p->ctrl_im.kp_flux = 2.0;
    p->ctrl_im.ki_flux = 20.0;
    p->ctrl_im.id_limit = 20.0;
    p->ctrl_im.iq_limit = 20.0;
    p->ctrl_im.vdq_limit = 520.0;
    p->ctrl_im.psi_r_ref = 0.75;
}
