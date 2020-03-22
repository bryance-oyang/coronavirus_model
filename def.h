/** @file
 * @brief parameter definitions
 */
#ifndef DEF_H
#define DEF_H

#include <stddef.h>

/* openmp parallelize parameters */
#define NTHREAD 4
#define THREAD_SCHEDULE static

/** total population */
#define POPULATION 4e8
/** number infected at start */
#define INFECTED_SEED 100

/* disease parameters */
/* assumes people are only infectious in a time window */
/** start of contagious phase (days after infected) */
#define CONTAG_T0 2
/** end of contagious phase, due to self-quarantine after recognizing
 * symptoms, etc (days after infected) */
#define CONTAG_T1 5
/** base reproductive number */
#define R0 2.2
/** reproductive number after mitigation */
#define MITIGATED_R0 0.8
/** time when mitigation takes effect (days) */
#define MITIGATION_START_TIME 25
/** fatality rate */
#define FATALITY_RATE 0.01
/** days until infection resolves (either recover or dead) */
#define NDAY_INFECTED 25

/** simulation time in days */
#define SIM_TIME 180

/* output */
#define MAX_OUTPUT 200

/** infected curve time resolution */
#define NCELL 1000

/** whether to do first order reconstruction */
#define RECONSTRUCT 1

/* timing */
/** actual dt is min of SUGGESTED_DT and numerically stable one */
#define SUGGESTED_DT 1000
/** Courant–Friedrichs–Lewy number */
#define CFL_NUM 0.36
/** max number of timesteps to take */
#define MAX_EPOCH 1e12

#define BUFLEN 1024

/* autocomputed */
#define OUTPUT_DT ((double) (SIM_TIME) / (MAX_OUTPUT))
#define NDAYPERCELL ((double)(NDAY_INFECTED) / (NCELL))

#endif /* DEF_H */
