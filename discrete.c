#include "discrete.h"
#include "def.h"

extern double dt;

/**
 * @brief compute d_vulnerable_dt
 */
void discrete_compute_vulnerable_slope(struct discrete *d, double
		ncontagious)
{
	double vulnerable_frac;

	vulnerable_frac = d->vulnerable / POPULATION;
	d->d_vulnerable_dt = -1 * vulnerable_frac * R0 * ncontagious
		/ (CONTAG_T1 - CONTAG_T0);
}

/**
 * @brief take timestep
 *
 * @param half_timestep 1 means this is half timestep, 0 means full
 * timestep
 */
void discrete_advance(struct discrete *d, int half_timestep)
{
	double step_dt;

	if (half_timestep) {
		step_dt = 0.5 * dt;
		d->vulnerable_copy = d->vulnerable;
		d->recovered_copy = d->recovered;
		d->dead_copy = d->dead;
	} else {
		step_dt = dt;
	}

	d->vulnerable = d->vulnerable_copy
		+ step_dt * d->d_vulnerable_dt;
	d->recovered = d->recovered_copy
		+ step_dt * d->d_recovered_dt;
	d->dead = d->dead_copy
		+ step_dt * d->d_dead_dt;
}
