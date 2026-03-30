#include "im_model.h"

#include <math.h>

#include "../control/transforms.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double dia;
    double dib;
    double dpsira;
    double dpsirb;
    double dtheta_m;
    double domega_m;
} ImDeriv;

static double wrap_pm_pi(double x)
{
    while (x > M_PI) x -= 2.0 * M_PI;
    while (x < -M_PI) x += 2.0 * M_PI;
    return x;
}

static ImDeriv im_deriv(const PlantState *x, const PlantInput *u, const ImMotorParams *mp)
{
    ImDeriv k;
    const double Tr = mp->Lr / mp->Rr;
    const double sigma = 1.0 - (mp->Lm * mp->Lm) / (mp->Ls * mp->Lr);
    const double wr = mp->pole_pairs * x->omega_m; /* rotor electrical speed */

    const double a = -(mp->Rs / (sigma * mp->Ls) + (mp->Lm * mp->Lm) / (sigma * mp->Ls * mp->Lr * Tr));
    const double b = 1.0 / (sigma * mp->Ls);
    const double c = mp->Lm / (sigma * mp->Ls * mp->Lr * Tr);
    const double d = mp->Lm / Tr;

    /* J*[x;y]=[-y;x] */
    const double J_psira = -x->psi_r_beta;
    const double J_psirb = x->psi_r_alpha;

    k.dia = a * x->i_alpha + c * x->psi_r_alpha - c * Tr * wr * J_psira + b * u->u_alpha;
    k.dib = a * x->i_beta  + c * x->psi_r_beta  - c * Tr * wr * J_psirb + b * u->u_beta;

    k.dpsira = d * x->i_alpha - (1.0 / Tr) * x->psi_r_alpha + wr * J_psira;
    k.dpsirb = d * x->i_beta  - (1.0 / Tr) * x->psi_r_beta  + wr * J_psirb;

    {
        const double Te = 1.5 * mp->pole_pairs * (mp->Lm / mp->Lr) *
                          (x->psi_r_alpha * x->i_beta - x->psi_r_beta * x->i_alpha);
        k.dtheta_m = x->omega_m;
        k.domega_m = (Te - u->T_load - mp->B * x->omega_m) / mp->J;
    }

    return k;
}

void im_init(PlantState *x)
{
    x->id = 0.0;
    x->iq = 0.0;
    x->i_alpha = 0.0;
    x->i_beta = 0.0;
    x->psi_r_alpha = 0.8;
    x->psi_r_beta = 0.0;
    x->theta_m = 0.0;
    x->omega_m = 0.0;
}

void im_step_rk4(PlantState *x, const PlantInput *u, const ImMotorParams *mp, double dt)
{
    ImDeriv k1 = im_deriv(x, u, mp);

    PlantState x2 = *x;
    x2.i_alpha += 0.5 * dt * k1.dia;
    x2.i_beta += 0.5 * dt * k1.dib;
    x2.psi_r_alpha += 0.5 * dt * k1.dpsira;
    x2.psi_r_beta += 0.5 * dt * k1.dpsirb;
    x2.theta_m += 0.5 * dt * k1.dtheta_m;
    x2.omega_m += 0.5 * dt * k1.domega_m;
    ImDeriv k2 = im_deriv(&x2, u, mp);

    PlantState x3 = *x;
    x3.i_alpha += 0.5 * dt * k2.dia;
    x3.i_beta += 0.5 * dt * k2.dib;
    x3.psi_r_alpha += 0.5 * dt * k2.dpsira;
    x3.psi_r_beta += 0.5 * dt * k2.dpsirb;
    x3.theta_m += 0.5 * dt * k2.dtheta_m;
    x3.omega_m += 0.5 * dt * k2.domega_m;
    ImDeriv k3 = im_deriv(&x3, u, mp);

    PlantState x4 = *x;
    x4.i_alpha += dt * k3.dia;
    x4.i_beta += dt * k3.dib;
    x4.psi_r_alpha += dt * k3.dpsira;
    x4.psi_r_beta += dt * k3.dpsirb;
    x4.theta_m += dt * k3.dtheta_m;
    x4.omega_m += dt * k3.domega_m;
    ImDeriv k4 = im_deriv(&x4, u, mp);

    x->i_alpha += dt * (k1.dia + 2.0 * k2.dia + 2.0 * k3.dia + k4.dia) / 6.0;
    x->i_beta += dt * (k1.dib + 2.0 * k2.dib + 2.0 * k3.dib + k4.dib) / 6.0;
    x->psi_r_alpha += dt * (k1.dpsira + 2.0 * k2.dpsira + 2.0 * k3.dpsira + k4.dpsira) / 6.0;
    x->psi_r_beta += dt * (k1.dpsirb + 2.0 * k2.dpsirb + 2.0 * k3.dpsirb + k4.dpsirb) / 6.0;
    x->theta_m += dt * (k1.dtheta_m + 2.0 * k2.dtheta_m + 2.0 * k3.dtheta_m + k4.dtheta_m) / 6.0;
    x->omega_m += dt * (k1.domega_m + 2.0 * k2.domega_m + 2.0 * k3.domega_m + k4.domega_m) / 6.0;

    x->theta_m = wrap_pm_pi(x->theta_m);
}

void im_get_output(const PlantState *x, const PlantInput *u, const ImMotorParams *mp, PlantOutput *y)
{
    double theta_flux = atan2(x->psi_r_beta, x->psi_r_alpha);
    double omega_e = mp->pole_pairs * x->omega_m;
    double id, iq;

    (void)u;

    park(x->i_alpha, x->i_beta, theta_flux, &id, &iq);

    y->i_alpha = x->i_alpha;
    y->i_beta = x->i_beta;
    inv_clarke(y->i_alpha, y->i_beta, &y->ia, &y->ib, &y->ic);

    y->id = id;
    y->iq = iq;

    y->psi_r_alpha = x->psi_r_alpha;
    y->psi_r_beta = x->psi_r_beta;
    y->psi_r_mag = hypot(x->psi_r_alpha, x->psi_r_beta);

    y->theta_m = x->theta_m;
    y->omega_m = x->omega_m;
    y->theta_e = theta_flux;
    y->omega_e = omega_e;

    y->Te = 1.5 * mp->pole_pairs * (mp->Lm / mp->Lr) *
            (x->psi_r_alpha * x->i_beta - x->psi_r_beta * x->i_alpha);

    {
        double eps = 1e-6;
        y->slip_e = (mp->Rr / (mp->Lr * (y->psi_r_mag + eps))) * iq * mp->Lm;
    }
}
