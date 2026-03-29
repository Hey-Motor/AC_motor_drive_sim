#ifndef ADC_MODEL_H
#define ADC_MODEL_H

#include "../core/types.h"

void adc_sample_currents(
    double ia_true,
    double ib_true,
    const AdcParams *ap,
    AdcSample *s
);

#endif