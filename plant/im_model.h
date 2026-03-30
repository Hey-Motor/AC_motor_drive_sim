#ifndef IM_MODEL_H
#define IM_MODEL_H

#include "../core/types.h"

void im_init(PlantState *x);
void im_step_rk4(PlantState *x, const PlantInput *u, const ImMotorParams *mp, double dt);
void im_get_output(const PlantState *x, const PlantInput *u, const ImMotorParams *mp, PlantOutput *y);

#endif
