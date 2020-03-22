#ifndef DEF_H
#define DEF_H

#include <stddef.h>

/* openmp parallelize parameters */
#define NTHREAD 4
#define THREAD_SCHEDULE static

/** total population */
#define POPULATION 4e8
/** number infected at start */
#define INFECTED_SEED 1000

/* disease parameters */
/** start of contagious phase after infected */
#define CONTAG_T0 2
/** end of contagious phase after infected */
#define CONTAG_T1 5
/** base reproductive number */
#define R0 2.2
/** fatality rate */
#define FATALITY_RATE 0.01

/** simulation time in days */
#define SIM_TIME 365

/* output */
#define MAX_OUTPUT 50

/* infected resolution */
#define NCELL 1000
#define NDAY 25

/** whether to do first order reconstruction */
#define RECONSTRUCT 1

/* timing */
#define SUGGESTED_DT 1000
#define CFL_NUM 0.36
#define MAX_EPOCH 1e12

#define BUFLEN 1024

/* autocomputed */
#define OUTPUT_DT ((double) (SIM_TIME) / (MAX_OUTPUT))
#define NDAYPERCELL ((double)(NDAY) / (NCELL))

#endif /* DEF_H */
