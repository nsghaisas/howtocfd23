
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>

#define f(x) sin(x)

void main(){

  FILE *fptr;

  int N = 1024,i;
  double pi,lx,dx,Nx=N;

  pi = 4*atan(1.0);
  lx = 2*pi;
  dx = lx/Nx;

  double *in_RHS, *out_final;
  fftw_complex *out_RHS,*in_RHS_modified;
  fftw_plan pdirect;
  fftw_plan pinverse;

  in_RHS = (double*) malloc(sizeof(double)*N);
  out_RHS = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  in_RHS_modified = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  out_final = (double*) malloc(sizeof(double)*N);

  pdirect = fftw_plan_dft_r2c_1d(N, in_RHS, out_RHS, FFTW_ESTIMATE);

  //Assign Values to the RHS

  for (i=0;i<N;i++)
    in_RHS[i] = f(i*dx);

  fftw_execute(pdirect);

  pinverse = fftw_plan_dft_c2r_1d(N, in_RHS_modified, out_final, FFTW_ESTIMATE);

  //Set Values of Fourier[f]/-s^2 except at x = 0,N/2 which is given by boundary conditions

  for (i=0;i<N/2;i++){
    if (i==0){
      in_RHS_modified[i][0] = 0;
      in_RHS_modified[i][1] = 0;
    }
    else{
      in_RHS_modified[i][0] = -out_RHS[i][0]/(i*i);
      in_RHS_modified[i][1] = -out_RHS[i][1]/(i*i);
    }
  }

  fftw_execute(pinverse);

  //OUTPUT

  fptr = fopen("poisson_1D_real_values.dat","w");
  //fprintf(fptr, "x f(x) phi(x)\n");
  for (i=0;i<N;i++)
      fprintf(fptr,"%g %g %g\n",i*dx,in_RHS[i],out_final[i]/Nx);
  fclose(fptr);

  /*printf("Input Values\n");
  for (i=0;i<N;i++)
      printf("%d - %f\n",i,in_RHS[i]);

  printf("\n\n");

  printf("Output Values\n");
  for (i=0;i<N;i++)
      printf("%d - %f\n",i,out_final[i]/Nx);
  printf("\n");*/


  fftw_destroy_plan(pdirect);
  fftw_destroy_plan(pinverse);
  fftw_free(in_RHS); fftw_free(out_RHS);
  fftw_free(in_RHS_modified);fftw_free(out_final);
  printf("DONE!!\nCheck Files for Output.\n");
}
