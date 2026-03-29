#ifndef PMSM_MODEL_H
#define PMSM_MODEL_H

#include "../core/types.h"

void pmsm_init(PlantState *x);
void pmsm_step_rk4(PlantState *x, const PlantInput *u, const MotorParams *mp, double dt);
void pmsm_get_output(const PlantState *x, const PlantInput *u, const MotorParams *mp, PlantOutput *y);

#endif