/****************************************************************
 *
 * diff_evo.c: Implementation of the differential evolution
 *	algorithm for global optimization
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
 *
 * https://www.potfit.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include "potfit.h"

#if defined(EVO)

#include "force.h"
#include "memory.h"
#include "optimize.h"
#include "potential_output.h"
#include "random.h"
#include "rescale.h"
#include "utils.h"

#define D (g_calc.ndimtot + 2)
#define NP 15 * D  // number of total population

#define JR 0.6  // jumping rate for opposite algorithm

// boundary values for self-adapting parameters
#define F_LOWER 0.1 /* lower value for F */
#define F_UPPER 0.9 /* upper value for F */
#define TAU_1 0.1   /* probability for changing F */
#define TAU_2 0.1   /* probability for changing CR */

/****************************************************************
 *
 *  initialize population with random numbers
 *
 ****************************************************************/

void init_population(double** pop, double* xi, double* cost) {
  // copy initial population into all populations
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < (D - 2); j++) pop[i][j] = xi[j];
    pop[i][D - 2] = F_LOWER + eqdist() * F_UPPER;
    pop[i][D - 1] = eqdist();
  }

  // create random populations (except for first one)
  for (int i = 1; i < NP; i++) {
    for (int j = 0; j < g_calc.ndim; j++) {
      double val = xi[g_pot.opt_pot.idx[j]];
      pop[i][g_pot.opt_pot.idx[j]] = 10.0 * val * (2 * eqdist() - 1);
    }
  }

  double forces[g_calc.mdim];

  for (int i = 0; i < NP; i++) cost[i] = calc_forces(pop[i], forces, 0);
}

/****************************************************************
 *
 *  differential evolution
 *
 ****************************************************************/

void run_differential_evolution(double* xi) {
  int a, b;                /* store randomly picked numbers */
                           //  int 	c; 			/* additional vector */
                           //  int 	d; 			/* additional vector */
                           //  int 	e; 			/* additional vector */
  int count = 0;           /* counter for loops */
  double cost_sum = 0.0;   /* average sum of squares for all configurations */
  double crit = 1000.0;    /* treshold for stopping criterion */
  double min_cost = 10e10; /* current minimum for all configurations */
  double max_cost = 0.0;   /* current maximum for all configurations */

  if (g_param.evo_threshold == 0.0) return;

  // vector for force calculation
  double* forces = (double*)Malloc(g_calc.mdim * sizeof(double));

  // vector with new configuration
  double* trial = (double*)Malloc(D * sizeof(double));
  memcpy(trial, xi, (D - 2) * sizeof(double));

  // allocate memory for all configurations
  double** pop_1 = (double**)Malloc(NP * sizeof(double*));
  double** pop_2 = (double**)Malloc(NP * sizeof(double*));
  double* best = (double*)Malloc(D * sizeof(double));
  double* cost = (double*)Malloc(NP * sizeof(double));

  for (int i = 0; i < NP; i++) {
    pop_1[i] = (double*)Malloc(D * sizeof(double));
    pop_2[i] = (double*)Malloc(D * sizeof(double));
  }

  init_population(pop_1, xi, cost);

  for (int i = 0; i < NP; i++) {
    if (cost[i] < min_cost) {
      min_cost = cost[i];
      memcpy(best, pop_1[i], D * sizeof(double));
    }
    if (cost[i] > max_cost) max_cost = cost[i];

    cost_sum += cost[i];
  }

  printf("done\n");

  crit = max_cost - min_cost;

  printf("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count, min_cost, cost_sum / (NP), crit);
  fflush(stdout);

  // main differential evolution loop
  while (crit >= g_param.evo_threshold && min_cost >= g_param.evo_threshold) {
    max_cost = 0.0;

    // randomly create new populations
    for (int i = 0; i < NP; i++) {
      // generate random numbers
      do
        a = (int)floor(eqdist() * NP);
      while (a == i);

      do
        b = (int)floor(eqdist() * NP);
      while (b == i || b == a);

      //       do
      //         c = (int)floor(eqdist() * NP);
      //       while (c == i || c == a || c == b);
      //
      //       do
      //         d = (int)floor(eqdist() * NP);
      //       while (d == i || d == a || d == b || d == c);
      //
      //       do
      //         e = (int)floor(eqdist() * NP);
      //       while (e == i || e == a || e == b || e == c || e == d);

      int j = (int)floor(eqdist() * g_calc.ndim);

      // self-adaptive parameters
      if (eqdist() < TAU_1)
        trial[D - 2] = F_LOWER + eqdist() * F_UPPER;
      else
        trial[D - 2] = pop_1[i][D - 2];

      if (eqdist() < TAU_2)
        trial[D - 1] = eqdist();
      else
        trial[D - 1] = pop_1[i][D - 1];

      double temp = 0.0;

      // create trail vectors with different methods
      for (int k = 1; k <= g_calc.ndim; k++) {
        if (eqdist() < trial[D - 1] || k == j) {
          /* DE/rand/1/exp */
          //           temp = pop_1[c][g_pot.opt_pot.idx[j]] + trial[D - 2] *
          //           (pop_1[a][g_pot.opt_pot.idx[j]] -
          //           pop_1[b][g_pot.opt_pot.idx[j]]);
          /* DE/best/1/exp */
          temp = best[g_pot.opt_pot.idx[j]] +
                 trial[D - 2] * (pop_1[a][g_pot.opt_pot.idx[j]] -
                                 pop_1[b][g_pot.opt_pot.idx[j]]);
          /* DE/rand/2/exp */
          //           temp = pop_1[e][j] + trial[D-2] * (pop_1[a][j] +
          //           pop_1[b][j] - pop_1[c][j] - pop_1[d][j]);
          /* DE/best/2/exp */
          //           temp = best[j] + trial[D-2] * (pop_1[a][j] + pop_1[b][j]
          //           - pop_1[c][j] - pop_1[d][j]);
          /* DE/rand-to-best/1/exp */
          //           temp = pop_1[c][j] + (1 - trial[D-2]) * (best[j] -
          //           pop_1[c][j]) + trial[D-2]
          //           * (pop_1[a][j] - pop_1[b][j]);
          /* DE/rand-to-best/2/exp */
          //           temp = pop_1[e][j] + (1 - trial[D-2]) * (best[j] -
          //           pop_1[e][j]) + trial[D-2]
          //           * (pop_1[a][j] + pop_1[b][j] - pop_1[c][j] -
          //           pop_1[d][j]);

          trial[g_pot.opt_pot.idx[j]] = temp;
        } else {
          trial[g_pot.opt_pot.idx[j]] = pop_1[i][g_pot.opt_pot.idx[j]];
        }

        j = (j + 1) % g_calc.ndim;
      }

      double force = calc_forces(trial, forces, 0);

      if (force < min_cost) {
        memcpy(best, trial, D * sizeof(double));

        if (*g_files.tempfile != '\0') {
          for (j = 0; j < g_calc.ndim; j++)

            xi[g_pot.opt_pot.idx[j]] = trial[g_pot.opt_pot.idx[j]];
          write_pot_table_potfit(g_files.tempfile);
        }
        min_cost = force;
      }

      if (force <= cost[i]) {
        memcpy(pop_2[i], trial, D * sizeof(double));

        cost[i] = force;

        if (force > max_cost) max_cost = force;
      } else {
        memcpy(pop_2[i], pop_1[i], D * sizeof(double));

        if (cost[i] > max_cost) max_cost = cost[i];
      }
    }

    cost_sum = 0.0;

    for (int i = 0; i < NP; i++) cost_sum += cost[i];

    printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count + 1, min_cost, cost_sum / (NP),
           max_cost - min_cost);
    fflush(stdout);

    for (int i = 0; i < NP; i++) memcpy(pop_1[i], pop_2[i], D * sizeof(double));
    count++;

    /* End optimization if break flagfile exists */
    if (*g_files.flagfile != '\0') {
      FILE* ff = fopen(g_files.flagfile, "r");

      if (NULL != ff) {
        printf("\nEvolutionary algorithm terminated ");
        printf("in presence of break flagfile \"%s\"!\n\n", g_files.flagfile);
        fclose(ff);
        remove(g_files.flagfile);
        break;
      }
    }

    crit = max_cost - min_cost;
  }

  printf("Finished differential evolution.\n");
  fflush(stdout);

  memcpy(xi, best, g_calc.ndimtot * sizeof(double));
}

#endif  // EVO
