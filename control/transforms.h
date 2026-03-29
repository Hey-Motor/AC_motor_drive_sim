#ifndef TRANSFORMS_H
#define TRANSFORMS_H

void clarke_2shunt(double ia, double ib, double *alpha, double *beta);
void inv_clarke(double alpha, double beta, double *ia, double *ib, double *ic);

void park(double alpha, double beta, double theta, double *d, double *q);
void inv_park(double d, double q, double theta, double *alpha, double *beta);

int circle_limit(double *x, double *y, double limit);

#endif