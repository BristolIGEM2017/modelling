#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<sys/time.h>
#include<omp.h>

int main(int argc, char *argv[]) {

  float dx, dt, final_t, a, v, diff, adv, err, tol;
  int npts, t;
  FILE *fp;

  // Simulation parameters
  npts = 101;
  final_t = 1000;
  dx = 1.;
  dt = 1.;
  tol = 1E-6;

  // Advection-Diffusion Coefficients
  v = 0.3;  // Transport velocity
  a = 0.15;  // Diffusion Coefficient

  // Set up of 1D array
  double *U, *Unew;
  U = (double *)calloc(npts, sizeof(double *) );
  Unew = (double *)calloc(npts, sizeof(double *) );

  // Check CFL condition
  float CFL = (2 * dt * a) / (dx * dx);
  if (CFL > 1.0) {
    printf("CFL of %f too large, reduce dt, diffusivity or increase dx\n", CFL);
    return EXIT_FAILURE;
  }

  // Impose initial conditions
  U[20] = 1.;
  U[30] = 2.;
  t = 0;

  // March in time until convergence
  err = 1E-4;
  while (err > tol) {

    // Set the convergence error to zero
    err = 0;

    // March in space
    #pragma omp parallel
    for (int i = 1; i < npts-1; i++) {
      diff = dt * a * ((U[i+1] + U[i-1] - 2 * U[i]) / (dx * dx));
      adv = dt * v * ((U[i+1] - U[i-1]) / (2 * dx));

      Unew[i] = U[i] + diff - adv;

      // Re-impose initial conditions in next timestep for error calculation
      Unew[20] = 1.;
      Unew[30] = 2.;
      if (err < fabs(Unew[i] - U[i])) {err = fabs(Unew[i] - U[i]);}
    };

    for (int i = 1; i < npts-1; i++) {
      U[i] = Unew[i];
    }

    t++;
  };

  // Data output
  fp = fopen("solution.dat", "w");

  fprintf(fp, "Final Time [s], Error, Values\n");
  fprintf(fp, "%15.8E, %15.8E" , t * dt, err);

  for (int i = 0; i < npts; i++) {
    fprintf(fp, ", %15.8E", U[i]);
  }

  fprintf(fp, "\n");

  // free memory
  fclose(fp);
  free(U);
  free(Unew);

  return EXIT_SUCCESS;
}
