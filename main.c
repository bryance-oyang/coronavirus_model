#include "def.h"
#include "infected.h"
#include "discrete.h"
#include <string.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>

double time;
double dt;

/* initial condition */
static void init(struct infected *inf, struct discrete *d)
{
	inf->n[0] = INFECTED_SEED / NDAYPERCELL;

	d->vulnerable = POPULATION - INFECTED_SEED;
	d->recovered = 0;
	d->dead = 0;

	infected_compute_total_infected(inf);
}

/**
 * @brief do second order time integration
 */
void advance_timestep(struct infected *inf, struct discrete *d)
{
	int half_timestep;
	double ncontagious;
	double d_resolve_dt;

	/* second order time integration */
	for (half_timestep = 1; half_timestep >= 0; half_timestep--) {
		/* new infections */
		ncontagious = infected_get_ncontagious(inf);
		discrete_compute_vulnerable_slope(d, ncontagious);

		/* infected */
		infected_set_left_boundary_J(inf, -1 * d->d_vulnerable_dt);
		infected_advance(inf, half_timestep);
		d_resolve_dt = infected_compute_resolve_rate(inf);

		/* discrete */
		d->d_recovered_dt = d_resolve_dt * (1 - FATALITY_RATE);
		d->d_dead_dt = d_resolve_dt * FATALITY_RATE;
		discrete_advance(d, half_timestep);
	}

	infected_compute_total_infected(inf);
}

/* output data to textfiles */
static void output(int epoch, struct infected *inf, struct discrete *d)
{
	int i;
	FILE *f;
	char filename[BUFLEN];

	/* output infected */
	sprintf(filename, "data/infected_%05d.dat", epoch);

	if ((f = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "fopen failed: %s\n", filename);
		raise(SIGSEGV);
	}

	for (i = 0; i < NCELL; i++) {
		fprintf(f, "%g %g\n", i*NDAYPERCELL, inf->n[i]);
	}

	fclose(f);

	/* output discrete */
	sprintf(filename, "data/discrete_%05d.dat", epoch);

	if ((f = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "fopen failed: %s\n", filename);
		raise(SIGSEGV);
	}

	fprintf(f, "%g\n", time);
	fprintf(f, "%g\n", (double)POPULATION + INFECTED_SEED);
	fprintf(f, "%g\n", (double)CONTAG_T0);
	fprintf(f, "%g\n", (double)CONTAG_T1);
	fprintf(f, "%g\n", inf->ninfected);
	fprintf(f, "%g\n", d->vulnerable);
	fprintf(f, "%g\n", d->recovered);
	fprintf(f, "%g\n", d->dead);

	fclose(f);
}


int main()
{
	int epoch, noutput;
	double output_time;
	struct infected inf;
	struct discrete d;

	dt = fmin(CFL_NUM * NDAYPERCELL, SUGGESTED_DT);

	/* initial condition */
	noutput = 0;
	time = 0;
	output_time = 0;
	init(&inf, &d);

	/* timestepping */
	for (epoch = 0; epoch < MAX_EPOCH; epoch++) {
		printf("t = %g | dt = %g\n", time, dt);
		if (time >= output_time) {
			output(noutput, &inf, &d);
			output_time = time + OUTPUT_DT;
			noutput++;
			if (noutput >= MAX_OUTPUT) {
				break;
			}
		}

		/* advance timestep */
		advance_timestep(&inf, &d);

		time += dt;
	}

	return 0;
}
