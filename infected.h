#ifndef INFECTED_H
#define INFECTED_H

#include "def.h"

struct infected {
	/** people per day */
	double n[NCELL];
	/** copy of n from start of timestep (for second order time
	 * integration) */
	double n0_copy[NCELL];
	/** current */
	double J[NCELL + 1];
	/** change per time from source term */
	double src[NCELL];
	/** reconstructed left state (indexed by boundary, left of boundary) */
	double nL[NCELL + 1];
	/** total number of people infected */
	double ninfected;
};

void infected_compute_total_infected(struct infected *inf);
void infected_advance(struct infected *inf, int half_timestep);
void infected_set_left_boundary_J(struct infected *inf, double J);
double infected_compute_resolve_rate(struct infected *inf);
double infected_get_ncontagious(struct infected *inf);

#endif /* INFECTED_H */
