#include "pmsm_model.h"
#include "../control/transforms.h"

typedef struct {
    double did;
    double diq;
    double dtheta_m;
    double domega_m;
} Deriv;

static Deriv pmsm_deriv(const PlantState *x, const PlantInput *u, const MotorParams *mp)
{
    Deriv k;
    double theta_e = mp->pole_pairs * x->theta_m;
    double omega_e = mp->pole_pairs * x->omega_m;

    double ud, uq;
    park(u->u_alpha, u->u_beta, theta_e, &ud, &uq);

    k.did = (ud - mp->Rs * x->id + omega_e * mp->Lq * x->iq) / mp->Ld;
    k.diq = (uq - mp->Rs * x->iq - omega_e * (mp->Ld * x->id + mp->psi_f)) / mp->Lq;

    double Te = 1.5 * mp->pole_pairs *
                (mp->psi_f * x->iq + (mp->Ld - mp->Lq) * x->id * x->iq);

    k.dtheta_m = x->omega_m;
    k.domega_m = (Te - u->T_load - mp->B * x->omega_m) / mp->J;
    return k;
}

void pmsm_init(PlantState *x)
{
    x->id = 0.0;
    x->iq = 0.0;
    x->theta_m = 3.14*20.0/180.0;
    x->omega_m = 0.0;
}

void pmsm_step_rk4(PlantState *x, const PlantInput *u, const MotorParams *mp, double dt)
{
    Deriv k1 = pmsm_deriv(x, u, mp);

    PlantState x2 = *x;
    x2.id      += 0.5 * dt * k1.did;
    x2.iq      += 0.5 * dt * k1.diq;
    x2.theta_m += 0.5 * dt * k1.dtheta_m;
    x2.omega_m += 0.5 * dt * k1.domega_m;
    Deriv k2 = pmsm_deriv(&x2, u, mp);

    PlantState x3 = *x;
    x3.id      += 0.5 * dt * k2.did;
    x3.iq      += 0.5 * dt * k2.diq;
    x3.theta_m += 0.5 * dt * k2.dtheta_m;
    x3.omega_m += 0.5 * dt * k2.domega_m;
    Deriv k3 = pmsm_deriv(&x3, u, mp);

    PlantState x4 = *x;
    x4.id      += dt * k3.did;
    x4.iq      += dt * k3.diq;
    x4.theta_m += dt * k3.dtheta_m;
    x4.omega_m += dt * k3.domega_m;
    Deriv k4 = pmsm_deriv(&x4, u, mp);

    x->id      += dt * (k1.did + 2.0 * k2.did + 2.0 * k3.did + k4.did) / 6.0;
    x->iq      += dt * (k1.diq + 2.0 * k2.diq + 2.0 * k3.diq + k4.diq) / 6.0;
    x->theta_m += dt * (k1.dtheta_m + 2.0 * k2.dtheta_m + 2.0 * k3.dtheta_m + k4.dtheta_m) / 6.0;
    x->omega_m += dt * (k1.domega_m + 2.0 * k2.domega_m + 2.0 * k3.domega_m + k4.domega_m) / 6.0;
}

void pmsm_get_output(const PlantState *x, const PlantInput *u, const MotorParams *mp, PlantOutput *y)
{
    (void)u;

    y->id = x->id;
    y->iq = x->iq;
    y->theta_m = x->theta_m;
    y->omega_m = x->omega_m;

    y->theta_e = mp->pole_pairs * x->theta_m;
    y->omega_e = mp->pole_pairs * x->omega_m;

    inv_park(x->id, x->iq, y->theta_e, &y->i_alpha, &y->i_beta);
    inv_clarke(y->i_alpha, y->i_beta, &y->ia, &y->ib, &y->ic);

    y->Te = 1.5 * mp->pole_pairs *
            (mp->psi_f * x->iq + (mp->Ld - mp->Lq) * x->id * x->iq);
}