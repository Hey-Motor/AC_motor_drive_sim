#include "adc_model.h"
#include "../control/transforms.h"
#include <math.h>

static double quantize(double x, double lsb)
{
    if (lsb <= 0.0) return x;
    return lsb * floor(x / lsb + 0.5);
}

void adc_sample_currents(
    double ia_true,
    double ib_true,
    const AdcParams *ap,
    AdcSample *s
)
{
    double ia = ia_true;
    double ib = ib_true;

    if (ap->enable_offset) {
        ia += ap->ia_offset;
        ib += ap->ib_offset;
    }

    if (ap->enable_quant) {
        ia = quantize(ia, ap->adc_lsb);
        ib = quantize(ib, ap->adc_lsb);
    }

    s->ia = ia;
    s->ib = ib;
    clarke_2shunt(ia, ib, &s->i_alpha, &s->i_beta);
}