#ifndef CONFIG_H
#define CONFIG_H

/* =========================
 * 仿真主配置
 * ========================= */
#define CTRL_FREQ_HZ        10000.0   /* 控制频率 10 kHz */
#define PLANT_SUBSTEPS      50        /* 每个控制周期内，连续对象积分子步数 */
#define LOG_DECIM           1         /* 每个控制周期都记录一次 */

/* =========================
 * 电机/对象模型类型
 * ========================= */
#define MOTOR_TYPE_PMSM     0
#define MOTOR_TYPE_IM       1
#define MOTOR_TYPE_DEFAULT  MOTOR_TYPE_IM

/* PMSM 子类型开关（仅在 MOTOR_TYPE_PMSM 时生效）
 * 1: IPMSM
 * 0: SPMSM
 */
#define USE_IPMSM           1

/* =========================
 * 逆变器模式
 * INV_IDEAL     : 理想电压直通，不做母线限制，不做非线性
 * INV_AVG_LINEAR: 平均模型 + 母线圆限幅
 * INV_AVG_NL    : 平均模型 + 母线圆限幅 + 非线性平均误差
 * ========================= */
#define INV_IDEAL           0
#define INV_AVG_LINEAR      1
#define INV_AVG_NL          2

/* 默认逆变器模式：后面你要做对比，改这里最方便 */
#define INVERTER_MODE_DEFAULT   INV_AVG_LINEAR

/* =========================
 * 采样与逆变器非理想开关
 * ========================= */
#define USE_ADC_MODEL       0
#define USE_INVERTER_NL     1

/* =========================
 * 观测器速度提取模式
 * OBS_SPEED_ATAN_DIFF : atan2 得角度，再微分+低通得速度
 * OBS_SPEED_PLL       : 磁链向量正交 PLL，同时输出角度和速度
 * ========================= */
#define OBS_SPEED_ATAN_DIFF     0
#define OBS_SPEED_PLL           1

/* 默认速度提取方式：推荐 PLL */
#define OBS_SPEED_MODE_DEFAULT  OBS_SPEED_PLL

/* =========================
 * 观测器模式选择
 * OBSERVER_FLUXPHI   : 你当前“非线性磁链+在线永磁磁链估计”那版
 * OBSERVER_PEBO_DREM : Bazylev/Pyrkin/Bobtsov 2018 这版
 * ========================= */
#define OBSERVER_FLUXPHI        0
#define OBSERVER_PEBO_DREM      1
#define OBSERVER_BOBTSOV2015    2
#define OBSERVER_IPMSM_NLRS     3
#define OBSERVER_INOUE2011      4
#define OBSERVER_BERNARD2021    5
#define OBSERVER_ETA_ONLY       6
#define OBSERVER_PIIPPO2008     7
/* 默认切到 Piippo2008 基础自适应观测器 */
#define OBSERVER_MODE_DEFAULT   OBSERVER_PIIPPO2008


#endif
