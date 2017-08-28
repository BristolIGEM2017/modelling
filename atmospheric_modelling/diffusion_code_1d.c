#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>
#include<omp.h>

int main(int argc, char *argv[]) {

  float dx, dt, final_t, a, v, diff, adv;
  int npts;
  FILE *fp;

  // Simulation parameters
  npts = 101;
  final_t = 1000;
  dx = 1;
  dt = 1;

  // Advection-Diffusion Coefficients
  v = 0.3;  // Transport velocity
  a = 0.15;  // Diffusion Coefficient

  // Set up of 1D array
  double *err;
  err = (double *)calloc((int)final_t/dt, sizeof(float *) );

  // Set up of 2D array
  double **U;
  U = (double **)calloc(npts, sizeof(double *) );

  for (int i = 0; i < npts; i++){
    U[i] = calloc((int)final_t/dt + 1, sizeof(double) );
  }

  // Check CFL condition
  float CFL = (2 * dt * a) / (dx * dx);
  if (CFL > 1.0) {
    printf("CFL of %f too large, reduce dt, diffusivity or increase dx\n", CFL);
    return EXIT_FAILURE;
  }

  // Impose initial conditions
  U[20][0] = 1.;
  U[30][0] = 2.;

  // March in time
  for (int t = 0; t < (int)(final_t/dt); t++) {

    // Set the convergence error to zero
    err[t] = 0;

    // March in space
    #pragma omp parallel
    for (int i = 1; i < npts-1; i++) {
      diff = dt * a * (U[i+1][t] + U[i-1][t] - 2 * U[i][t]) / (dx * dx);
      adv = dt * v * (U[i+1][t] - U[i-1][t]) / (2 * dx);

      U[i][t+1] = U[i][t] + diff - adv;

      // Re-impose initial conditions in next timestep for error calculation
      U[20][t+1] = 1.;
      U[30][t+1] = 2.;
      if (err[t] < fabs(U[i][t+1] - U[i][t])) {err[t] = fabs(U[i][t+1] - U[i][t]);}
    };
  };

  // Data output
  fp = fopen("solution.dat", "w");

  fprintf(fp, "Time [s], Max Error, Values\n");
  for (int t = 0; t < (int)(final_t/dt); t++) {
    fprintf(fp, "%15.8E, %15.8E" , t * dt, err[t]);
    for (int i = 0; i < npts; i++) {
      fprintf(fp, ", %15.8E", U[i][t]);
    }
    fprintf(fp, "\n");
  }

  // free memory
  fclose(fp);
  for(int i = 0; i < npts; i++) free(U[i]);

  return EXIT_SUCCESS;
}
