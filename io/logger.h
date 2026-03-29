#ifndef LOGGER_H
#define LOGGER_H

#include "../core/types.h"
#include <stdio.h>

typedef struct {
    FILE *fp;
} CsvLogger;

int logger_open(CsvLogger *lg, const char *filename);

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
);

void logger_close(CsvLogger *lg);

#endif