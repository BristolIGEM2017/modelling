#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>
#include<omp.h>

int main(int argc, char *argv[]) {

  float dx, dt, final_t, a, v, err, diff, adv;
  int npts;
  FILE *fp;

  // Simulation parameters
  npts = 101;
  final_t = 100;
  dx = 1;
  dt = 1;
  err = 0;

  // Advection-Diffusion Coefficients
  v = 0.3;  // Transport velocity
  a = 0.15;  // Diffusion Coefficient

  // Set up of 2D array
  double **U;
  U = (double **)calloc(npts, sizeof(double *) );

  for (int i = 0; i < npts; i++){
    U[i] = calloc((int)(final_t/dt), sizeof(double) );
  }

  // March in time
  for (int t = 0; t < (int)(final_t/dt); t++) {
    printf("Timestep: %d\n", t);

    // Impose initial conditions
    U[20][t] = 1.;
    U[30][t] = 2.;

    // March in space
    #pragma omp parallel
    for (int i = 1; i < npts-1; i++) {
      diff = dt * a * (U[i+1][t] + U[i-1][t] - 2 * U[i][t]) / (dx * dx);
      adv = dt * v * (U[i+1][t] - U[i-1][t]) / (2 * dx);

      U[i][t+1] = U[i][t] + diff - adv;
    };
  };

  // Data output
  fp = fopen("solution.dat", "w");

  fprintf(fp, "Time [s], Values\n");
  for (int t = 0; t < (int)(final_t/dt); t++) {
    for (int i = 0; i < npts; i++) {
      fprintf(fp, "%15.8E ,", U[i][t]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  free(U);

  return EXIT_SUCCESS;
}
