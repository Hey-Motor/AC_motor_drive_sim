#ifndef TYPES_H
#define TYPES_H

typedef struct {
    double Ts_ctrl;
    double dt_plant;
    double t_end;
    int substeps;
    int log_decim;
} SimConfig;

typedef struct {
    int pole_pairs;
    double Rs;
    double Ld;
    double Lq;
    double psi_f;
    double J;
    double B;
} MotorParams;

typedef struct {
    int pole_pairs;
    double Rs;
    double Rr;
    double Lm;
    double Ls;
    double Lr;
    double J;
    double B;
} ImMotorParams;

typedef struct {
    double Vdc;

    /* 非线性平均误差的饱和值，相对母线电压比例 */
    double nl_sat_ratio;

    /* 电流达到该阈值后，非线性误差进入饱和 */
    double nl_current_thr;

    /* 当前逆变器模式：INV_IDEAL / INV_AVG_LINEAR / INV_AVG_NL */
    int mode;

    /* 是否启用非线性 */
    int enable_nl;
} InverterParams;

typedef struct {
    double ia_offset;
    double ib_offset;
    double adc_lsb;
    int enable_offset;
    int enable_quant;
} AdcParams;

typedef struct {
    double kp_id, ki_id;
    double kp_iq, ki_iq;
    double kp_spd, ki_spd;
    double iq_limit;
    double vdq_limit;
} ControlParams;

typedef struct {
    double kp_id, ki_id;
    double kp_iq, ki_iq;
    double kp_spd, ki_spd;
    double kp_flux, ki_flux;
    double id_limit;
    double iq_limit;
    double vdq_limit;
    double psi_r_ref;
} ImControlParams;

typedef struct {
    SimConfig sim;
    MotorParams motor;
    ImMotorParams motor_im;
    InverterParams inverter;
    AdcParams adc;
    ControlParams ctrl;
    ImControlParams ctrl_im;
} Params;

/* =========================
 * 连续对象状态
 * 用 dq 状态积分
 * ========================= */
typedef struct {
    double id;
    double iq;
    double i_alpha;
    double i_beta;
    double psi_r_alpha;
    double psi_r_beta;
    double theta_m;
    double omega_m;
} PlantState;

/* =========================
 * 连续对象输入
 * 平均逆变器输出给电机的 alpha-beta 电压
 * ========================= */
typedef struct {
    double u_alpha;
    double u_beta;
    double T_load;
} PlantInput;

/* =========================
 * 连续对象输出
 * 既保留 dq 量，也保留三相量和 alpha-beta 量
 * 方便控制、观测器、画图统一使用
 * ========================= */
typedef struct {
    double ia, ib, ic;
    double i_alpha, i_beta;
    double id, iq;
    double theta_m, theta_e;
    double omega_m, omega_e;
    double Te;
    double psi_r_alpha, psi_r_beta;
    double psi_r_mag;
    double slip_e;
} PlantOutput;

typedef struct {
    double ia;
    double ib;
    double i_alpha;
    double i_beta;
} AdcSample;

/* =========================
 * 观测器输出
 * theta_e_hat / omega_e_hat : 给 FOC 用
 * 其余量保留出来，后面你做日志和调试会很方便
 * ========================= */
typedef struct {
    double theta_e_hat;
    double omega_e_hat;

    double psi_alpha_hat;
    double psi_beta_hat;
    double phi_hat;

    /* 新增：PEBO+DREM 调试量 */
    double R_hat;
    double eta1_hat;
    double eta2_hat;
    double beta_hat;
    double detQ;
        /* 原始 6 维回归 + 标量输出，专门用于诊断 DREM */
    double q1;
    double q2;
    double q3;
    double q4;
    double q5;
    double q6;
    double yreg;

        /* 原始中间量：用于检查论文式(11) */
    double z21;
    double z22;
    double xi1;
    double xi2;
} ObserverOutput;

/* =========================
 * 观测器内部状态
 * psi_alpha_hat, psi_beta_hat : 总磁链估计
 * phi_hat                     : 永磁磁链幅值估计
 *
 * theta_prev                  : atan2+微分法时用
 * omega_e_hat_f               : 速度滤波值
 *
 * pll_theta                   : PLL 内部角度
 * pll_omega_i                 : PLL 积分器状态
 *
 * initialized                 : 初始化标志
 * ========================= */
typedef struct {
    double psi_alpha_hat;
    double psi_beta_hat;
    double phi_hat;

    /* atan2 + 微分法所需状态 */
    double theta_prev;
    double omega_e_hat_f;

    /* PLL 所需状态 */
    double pll_theta;
    double pll_omega_i;

    int initialized;
} ObserverState;

typedef struct {
    double id_int;
    double iq_int;
    double spd_int;
    double flux_int;
} FocState;

#endif
