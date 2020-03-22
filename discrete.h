#ifndef DISCRETE_H
#define DISCRETE_H

struct discrete {
	double vulnerable;
	double recovered;
	double dead;

	double vulnerable_copy;
	double recovered_copy;
	double dead_copy;

	double d_vulnerable_dt;
	double d_recovered_dt;
	double d_dead_dt;
};

void discrete_compute_vulnerable_slope(struct discrete *d, double
		ncontagious);
void discrete_advance(struct discrete *d, int half_timestep);

#endif /* DISCRETE_H */
