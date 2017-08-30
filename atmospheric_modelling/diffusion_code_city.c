#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>

typedef struct {
    float *x;
    float *y;
    float *val;
} Source;

int count_sources(char file_name[]) {
  FILE* stream = fopen(file_name, "r");
  char line[64];
  int cnt = 0;
  if (stream != NULL) {
    while (fgets(line, 64, stream)) {
      cnt++;
    }
  } else {
    printf("File not found.");
  }
  fclose(stream);
  return cnt;
}

float read_sources(char file_name[], Source sources) {
  FILE* stream = fopen(file_name, "r");
  char line[64];
  int i = 0;
  if (stream != NULL) {
    while (fgets(line, 64, stream)) {
      fscanf(stream, "%f, %f, %f", &sources.x[i], &sources.y[i] ,&sources.val[i]);
      i++;
    }
  } else {
    printf("File not found.");
  }
  fclose(stream);
}

int main(int argc, char *argv[]) {

  float dx, dy, dt, final_t, a, u, v, diff, adv, err, tol;
  int nptsi, nptsj, t;
  FILE *fp;

  // Simulation parameters
  nptsi = 101;
  nptsj = 101;
  final_t = 1000;
  dx = 1.;
  dy = 1.;
  dt = 0.25;
  tol = 1E-6;

  // Read in sources file
  char source_file[] = "sources.dat";
  int lines = count_sources("sources.dat");

  // Create struct that will contain all sources and allocate memory
  Source sources;
  sources.x = calloc(lines, sizeof(double) );
  sources.y = calloc(lines, sizeof(double) );
  sources.val = calloc(lines, sizeof(double) );

  read_sources(source_file, sources);

  // Advection-Diffusion Coefficients
  u = 0.2;  // Transport velocity x
  v = 0.2;  // Transport velocity y
  a = 0.8;  // Diffusion Coefficient

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
  t = 0;

  // March in time until convergence
  err = 1E-4;
  while (err > tol) {

    // Set the convergence error to zero
    err = 0;

    // March in space
    #pragma omp parallel
    for (int i = 1; i < nptsi-1; i++) {
      for (int j = 1; j < nptsj-1; j++) {
        diff = (dt * a * (U[i+1][j] + U[i-1][j] - (2. * U[i][j])) / (dx * dx)) +
               (dt * a * (U[i][j+1] + U[i][j-1] - (2. * U[i][j])) / (dy * dy));
        adv = (dt * u * (U[i+1][j] - U[i-1][j]) / (2. * dx)) +
              (dt * v * (U[i][j+1] - U[i][j-1]) / (2. * dy));

        Unew[i][j] = U[i][j] + diff - adv;

        // Re-impose initial conditions
        Unew[25][25] = 2.;
        if (err < fabs(Unew[i][j] - U[i][j])) {err = fabs(Unew[i][j] - U[i][j]);}
      };
    };

    for (int i = 1; i < nptsi-1; i++) {
      for (int j = 1; j < nptsj-1; j++) {
        U[i][j] = Unew[i][j];
      };
    };

    t++;
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
  for(int i = 0; i < nptsi; i++) free(Unew[i]);

  return EXIT_SUCCESS;
}
