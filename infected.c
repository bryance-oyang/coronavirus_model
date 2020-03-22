#include "infected.h"
#include <math.h>

extern double dt;

/* for reconstruction */
static inline double slope_lim(double r)
{
	return fmax(0, fmin(1, r));
}

static void reconstruct(struct infected *inf)
{
	int i;

	/* loop over non-boundary cells */
#pragma omp parallel for simd schedule(THREAD_SCHEDULE) num_threads(NUM_THREAD)
	for (i = 1; i < NCELL-1; i++) {
		double n0, n1, n2, r;

		if (RECONSTRUCT) {
			n0 = inf->n[i-1];
			n1 = inf->n[i];
			n2 = inf->n[i+1];

			r = (n1 - n0) / (n2 - n1);

			inf->nL[i+1] = n1 + slope_lim(r) * (n2 - n1);
		} else {
			inf->nL[i+1] = inf->n[i];
		}
	}

	/* boundary */
	inf->nL[1] = inf->n[0];
	inf->nL[NCELL] = inf->n[NCELL-1];
}

static void calculate_J(struct infected *inf)
{
	int i;
	double speed;

	speed = 1;

	/* loop over faces */
#pragma omp simd
	for (i = 1; i < NCELL+1; i++) {
		inf->J[i] = inf->nL[i] * speed;
	}
}

static void calculate_src(struct infected *inf)
{
	int i;

#pragma omp simd
	for (i = 0; i < NCELL; i++) {
		inf->src[i] = 0;
	}
}

void infected_compute_total_infected(struct infected *inf)
{
	int i;

	/* calculate total infected */
	inf->ninfected = 0;
	for (i = 0; i < NCELL; i++) {
		inf->ninfected += inf->n[i] * NDAYPERCELL;
	}
}

/**
 * @brief take timestep
 *
 * @param half_timestep 1 means this is half timestep, 0 means full
 * timestep
 */
void infected_advance(struct infected *inf, int half_timestep)
{
	int i;
	double step_dt;
	double step_dt_dx;

	if (half_timestep) {
		/* determine dt */
		step_dt = 0.5 * dt;
		/* copy state */
#pragma omp simd
		for (i = 0; i < NCELL; i++) {
			inf->n0_copy[i] = inf->n[i];
		}
	} else {
		step_dt = dt;
	}

	/* calculate time derivatives */
	reconstruct(inf);
	calculate_J(inf);
	calculate_src(inf);

	/* add fluxes and src */
	step_dt_dx = step_dt  / NDAYPERCELL;
#pragma omp parallel for simd schedule(THREAD_SCHEDULE) num_threads(NUM_THREAD)
	for (i = 0; i < NCELL; i++) {
		inf->n[i] = inf->n0_copy[i] + step_dt_dx * (inf->J[i] - inf->J[i+1])
			+ step_dt * inf->src[i];
	}
}

void infected_set_left_boundary_J(struct infected *inf, double J)
{
	inf->J[0] = J;
}

/**
 * @brief compute number resolve per day (recover + dead)
 */
double infected_compute_resolve_rate(struct infected *inf)
{
	return inf->J[NCELL];
}

double infected_get_ncontagious(struct infected *inf)
{
	int i;
	double x;
	double N;

	N = 0;
	for (i = 0; i < NCELL; i++) {
		x = i * NDAYPERCELL;
		if (x >= CONTAG_T0 && x <= CONTAG_T1) {
			N += NDAYPERCELL * inf->n[i];
		}
	}

	return N;
}
