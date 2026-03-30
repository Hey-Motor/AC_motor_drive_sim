#include "logger.h"

int logger_open(CsvLogger *lg, const char *filename)
{
    lg->fp = fopen(filename, "w");
    if (!lg->fp) return 0;

    fprintf(lg->fp,
            "t,speed_ref,id_ref,iq_ref,"
            "u_alpha_cmd,u_beta_cmd,u_alpha_act,u_beta_act,"
            "ia,ib,ic,i_alpha,i_beta,id,iq,"
            "adc_ia,adc_ib,adc_i_alpha,adc_i_beta,"
            "theta_m,theta_e,omega_m,omega_e,"
            "psi_r_alpha,psi_r_beta,psi_r_mag,slip_e,"
            "theta_e_hat,omega_e_hat,theta_err,"
            "psi_alpha_hat,psi_beta_hat,phi_hat,"
            "R_hat,eta1_hat,eta2_hat,beta_hat,detQ,"
            "q1,q2,q3,q4,q5,q6,yreg,"
            "z21,z22,xi1,xi2,"
            "Te,T_load\n");
    return 1;
}

void logger_write(
    CsvLogger *lg,
    double t,
    double speed_ref,
    double id_ref,
    double iq_ref,
    double u_alpha_cmd,
    double u_beta_cmd,
    double u_alpha_act,
    double u_beta_act,
    const PlantOutput *y,
    const AdcSample *adc,
    const ObserverOutput *obs,
    double theta_err,
    double T_load
)
{
    if (!lg->fp) return;

    fprintf(lg->fp,
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f,%.9f,%.9f,"
        "%.9f,%.9f\n",
        t, speed_ref, id_ref, iq_ref,
        u_alpha_cmd, u_beta_cmd, u_alpha_act, u_beta_act,
        y->ia, y->ib, y->ic, y->i_alpha, y->i_beta, y->id, y->iq,
        adc->ia, adc->ib, adc->i_alpha, adc->i_beta,
        y->theta_m, y->theta_e, y->omega_m, y->omega_e,
        y->psi_r_alpha, y->psi_r_beta, y->psi_r_mag, y->slip_e,
        obs->theta_e_hat, obs->omega_e_hat, theta_err,
        obs->psi_alpha_hat, obs->psi_beta_hat, obs->phi_hat,
        obs->R_hat, obs->eta1_hat, obs->eta2_hat, obs->beta_hat, obs->detQ,
        obs->q1, obs->q2, obs->q3, obs->q4, obs->q5, obs->q6, obs->yreg,
        obs->z21, obs->z22, obs->xi1, obs->xi2,
        y->Te, T_load);
}

void logger_close(CsvLogger *lg)
{
    if (lg->fp) {
        fclose(lg->fp);
        lg->fp = NULL;
    }
}
