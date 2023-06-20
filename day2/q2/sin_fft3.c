#include<stdio.h>
#include <stdlib.h>
#include<math.h>
# include "fftw3.h"

void grid(int nx, double *x, double dx)
{
  int i;

  // populate x[i] s
  for(i=0; i<nx; i++)
    x[i] = ((double)i)  * dx;

  //// debug -- print x
  //printf("--x--\n");
  //for(i=0; i<nx; i++)
  //  printf("%d %lf\n", i, x[i]);
}

void setup_wavenos(int nx, double twopi_by_L, double *kx)
{
  int i;
  
  // setup the wavenumbers

  // first half are in regular order
  for(i=0; i<nx/2; i++)
    kx[i] = (double)i * twopi_by_L;

  // second half are in reverse order
  for(i=nx/2; i<nx; i++)
    kx[i] = (double)(-nx+i) * twopi_by_L;
}

void set_initial_fields(int nx, int ny, double *x, double *y, double **psi, double **omg)
{
  int i, j;

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      psi[i][j] = 0.0; // dummy init

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      omg[i][j] = sin(5.0*x[i]);
}

void output_soln(int nx, int ny, int iter, double *x, double *y, double **psi, double **omg)
{
  int i, j;
  FILE* fp;
  char fname[100];

  sprintf(fname, "psi_omg_%03d_%03d_%04d.dat", nx, ny, iter);

  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     fprintf(fp, "%lf %lf %lf %lf\n", x[i], y[j], psi[i][j], omg[i][j]);
  fclose(fp);

  printf("Done writing solution for stamp = %d to file %s\n", iter, fname);
}

int main()
{
  
  int nx, ny, i, ncx, ncy;
  double *x, dx, *y, dy, Lx, Ly, *kx, *ky, *f;
  double **psi, **omg;
  fftw_complex *fhat;
  fftw_plan plan_forward, plan_backward;

  nx = 128;      ny = nx;
  Lx = 2.0*M_PI; Ly = Lx;

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  dx = Lx/(double)nx;  dy = dx;
  x  = (double *)malloc(nx     * sizeof(double));   // grid points in x
  y  = (double *)malloc(ny     * sizeof(double));   // grid points in y
  kx = (double *)malloc(nx     * sizeof(double));   // wavenos in x
  ky = (double *)malloc(ny     * sizeof(double));   // wavenos in y
  printf("   >> Done allocating 1D arrays -- \n");

  // allocate 2D arrays dynamically
  // -- for streamfunction --
  psi = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    psi[i] = (double *)malloc(ny*sizeof(double));

  // -- for vorticity --
  omg = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    omg[i] = (double *)malloc(ny*sizeof(double));

  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");
 
  // define size of the transformed fields
  ncx = nx;        // first dimension is unchanged
  ncy = ny/2 + 1;  // second dimension is close to halved

  // allocate f and fhat using fftw_malloc 
  f = fftw_malloc(nx*ny*sizeof(double));
  fhat = fftw_malloc(ncx*(ncy+1)*sizeof(fftw_complex));

  // initialize the grid
  grid(nx, x, dx);  // -- along x --
  grid(ny, y, dy);  // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  setup_wavenos(nx, 1.0, kx);
  setup_wavenos(ny, 1.0, ky);

  set_initial_fields(nx,ny,x,y,omg,psi);  // initial fields
  output_soln(nx, ny, 0, x, y, psi, omg); // output initial omega fields

  // create plan for forward  Fourier transform
  plan_forward = fftw_plan_dft_r2c_2d(nx, ny, f, fhat, FFTW_ESTIMATE);

  // create plan for backward Fourier transform
  plan_backward = fftw_plan_dft_c2r_2d(nx, ny, fhat, f, FFTW_ESTIMATE);

  // execute forward Fourier transform of omega
    // pack omg into f
    fftw_execute(plan_forward);  // execute fftw
    // unpack fhat to omghat

  // solve for psihat

  // execute backward Fourier transform of psi and scale
    // pack psihat into fhat
    fftw_execute(plan_backward);  // execute fftw
    // unpack f to psi

  output_soln(nx, ny, 1, x, y, psi, omg); // output fields after operations

  // free memory
    fftw_free(fhat);     fftw_free(f);

   // ----1D arrays ---
   free(ky);  free(kx);  free(y);  free(x);
   // --- Done 1D arrays ---

   // ----2D arrays ---
   for(i=0; i<nx; i++)     free(omg[i]);    free(omg);
   for(i=0; i<nx; i++)     free(psi[i]);    free(psi);
  printf("\n > Done freeing up memory --------- \n");

  return 0;
}
