
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>


double f(double x,double y){
  return (sin(x)+sin(y));
}

void main(){

  FILE *fptr;

  int N = 256,i,j,Nh = N/2+1;
  double pi,lx,dx,Nx=N,ly,dy,Ny=N,k;

  pi = 4*atan(1.0);
  lx = 2*pi;  ly = 2*pi;
  dx = lx/Nx; dy = ly/Ny;

  double *omega, *psi;
  fftw_complex *ft_omega,*modded_ft_omega;
  fftw_plan pdirect;
  fftw_plan pinverse;

  omega = (double*) malloc(sizeof(double)*N*N);
  ft_omega = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N/2+1)*N);
  modded_ft_omega = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N/2+1)*N);
  psi = (double*) malloc(sizeof(double)*N*N);

  pdirect = fftw_plan_dft_r2c_2d(N, N, omega, ft_omega, FFTW_ESTIMATE);
  pinverse = fftw_plan_dft_c2r_2d(N, N, modded_ft_omega, psi, FFTW_ESTIMATE);

  //Assign Values to the RHS
  for (i=0;i<N;i++)
    for(j=0;j<N;j++)
      omega[i*N+j] = f(i*dx,j*dy);

  fftw_execute(pdirect);

  //Set Values of Fourier[f]/-s^2 except at s=0 which is given by boundary conditions
  for (i=0;i<N;i++){
    if(2*i<N)
      k = i*i;
    else
      k = (i-N)*(i-N);
    for (j=0;j<Nh;j++){
      double fac = -1.0*(k+j*j);
      if (fabs(fac)==0){
        modded_ft_omega[i*Nh+j][0] = 0;
        modded_ft_omega[i*Nh+j][1] = 0;
      }
      else{
        modded_ft_omega[i*Nh+j][0] = ft_omega[i*Nh+j][0]/fac;
        modded_ft_omega[i*Nh+j][1] = ft_omega[i*Nh+j][1]/fac;
      }
    }
  }

  fftw_execute(pinverse);

  //OUTPUT
  fptr = fopen("poisson_2D_real_values_c.dat","w");
  for (i=0;i<N;i++)
    for(j=0;j<N;j++)
      fprintf(fptr,"%lf %lf %lf %lf\n",i*dx,j*dy,omega[i*N+j],psi[i*N+j]/(Nx*Ny));
  fclose(fptr);


  fftw_destroy_plan(pdirect);
  fftw_destroy_plan(pinverse);
  fftw_free(omega); fftw_free(ft_omega);
  fftw_free(modded_ft_omega);fftw_free(psi);
  printf("\nDONE!!\nCheck Files for Output.\n");
}
