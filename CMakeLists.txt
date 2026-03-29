cmake_minimum_required(VERSION 3.15)
project(pmsm_obs_sim C)

set(CMAKE_C_STANDARD 99)

add_executable(pmsm_obs_sim
    main.c
    config/params.c
    plant/pmsm_model.c
    plant/inverter_avg.c
    control/transforms.c
    control/foc_controller.c
    control/observer_select.c
    control/observer_fluxphi.c
    control/observer_pebo_drem.c
    control/observer_bobtsov2015.c
    control/observer_ipmsm_nlrs.c
    io/adc_model.c
    control/observer_inoue2011.c
    control/observer_bernard2021.c
    control/observer_eta_only.c
    control/observer_piippo2008.c
    io/logger.c
)

target_link_libraries(pmsm_obs_sim m)
