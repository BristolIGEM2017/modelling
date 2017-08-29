#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>
#include<omp.h>

int main(int argc, char *argv[]) {

  float dx, dy, dt, final_t, a, v, u, diff, adv;
  int nptsi, nptsj;
  FILE *fp;

  // Simulation parameters
  nptsi = 101;
  nptsj = 101;
  final_t = 1000;
  dx = 1.;
  dt = 0.1;
  dy = 1.;

  // Advection-Diffusion Coefficients
  v = 0.2;  // Transport velocity x
  u = 0.2;  // Transport velocity y
  a = 0.8;  // Diffusion Coefficient

  // Set up of 1D array
  double *err;
  err = (double *)calloc((int)final_t/dt, sizeof(float *) );

  // Set up of 2D array
  double **U, **Unew;
  U = (double **)calloc(nptsi, sizeof(double *) );
  Unew = (double **)calloc(nptsi, sizeof(double *) );

  for (int i = 0; i < nptsj; i++){
    U[i] = calloc(nptsj, sizeof(double) );
    Unew[i] = calloc(nptsj, sizeof(double) );
  }

  // Check CFL condition
  float CFLx = (2 * dt * a) / (dx * dx);
  float CFLy = (2 * dt * a) / (dy * dy);
  if ((CFLx > 1.0) || (CFLy > 1.0)) {
    printf("CFL (%f/%f) too large, reduce dt, diffusivity or increase (dx/dy)\n", CFLx, CFLy);
    return EXIT_FAILURE;
  }

  // Impose initial conditions
  U[25][25] = 2.;

  // March in time
  for (int t = 0; t < (int)(final_t/dt); t++) {

    // Set the convergence error to zero
    err[t] = 0;

    // March in space
    #pragma omp parallel
    for (int i = 1; i < nptsi-1; i++) {
      for (int j = 1; j < nptsj-1; j++) {
        diff = (dt * a * (U[i+1][j] + U[i-1][j] - (2. * U[i][j])) / (dx * dx)) +
               (dt * a * (U[i][j+1] + U[i][j-1] - (2. * U[i][j])) / (dy * dy));
        adv = (dt * u * (U[i+1][j] - U[i-1][j]) / (2. * dx)) +
              (dt * v * (U[i][j+1] - U[i][j-1]) / (2. * dy));

        Unew[i][j] = U[i][j] + diff - adv;
      };
    };

    U = Unew;

    // Re-impose initial conditions
    U[25][25] = 2.;
  };


  // Data output
  fp = fopen("solution.dat", "w");

  for (int i = 0; i < nptsi; i++) {
    for (int j = 0; j < nptsj; j++) {
      fprintf(fp, "%15.8E,", U[i][j]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");


  // free memory
  fclose(fp);
  for(int i = 0; i < nptsi; i++) free(U[i]);

  return EXIT_SUCCESS;
}
