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
}