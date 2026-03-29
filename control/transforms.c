#include "transforms.h"
#include <math.h>

void clarke_2shunt(double ia, double ib, double *alpha, double *beta)
{
    const double inv_sqrt3 = 0.5773502691896258;
    *alpha = ia;
    *beta  = (ia + 2.0 * ib) * inv_sqrt3;
}

void inv_clarke(double alpha, double beta, double *ia, double *ib, double *ic)
{
    *ia = alpha;
    *ib = -0.5 * alpha + 0.8660254037844386 * beta;
    *ic = -0.5 * alpha - 0.8660254037844386 * beta;
}

void park(double alpha, double beta, double theta, double *d, double *q)
{
    double c = cos(theta);
    double s = sin(theta);
    *d =  c * alpha + s * beta;
    *q = -s * alpha + c * beta;
}

void inv_park(double d, double q, double theta, double *alpha, double *beta)
{
    double c = cos(theta);
    double s = sin(theta);
    *alpha = c * d - s * q;
    *beta  = s * d + c * q;
}

int circle_limit(double *x, double *y, double limit)
{
    double mag = sqrt((*x) * (*x) + (*y) * (*y));
    if (mag > limit && mag > 1e-12) {
        double k = limit / mag;
        *x *= k;
        *y *= k;
        return 1;
    }
    return 0;
}