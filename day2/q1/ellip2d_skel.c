#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//const double PI = 4.0*atan(1.0);

void grid(int nx, double xst, double xen, double *x, double *dx)
{
  int i;
  double dxunif;

  // uniform mesh for now; 
  // can use stretching factors to place x nodes later
  dxunif = (xen-xst)/(double)(nx-1);

  // populate x[i] s
  for(i=0; i<nx; i++)
    x[i] = xst + ((double)i)  * dxunif;

  // dx[i] s are spacing between adjacent xs
  for(i=0; i<nx-1; i++)
    dx[i] = ...;

  //// debug -- print x
  printf("--x--\n");
  for(i=0; i<nx; i++)
    printf("%d %lf\n", i, x[i]);

  //// debug -- print dx
  printf("--dx--\n");
  for(i=0; i<nx-1; i++)
    printf("%d %lf\n", i, dx[i]);
}

void set_initial_guess(int nx, int ny, double *x, double *y, double **T)
{
  int i, j;

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      T[i][j] = ...;
}

void calc_diffusivity(int nx, int ny, double *x, double *y, double **T, double **kdiff)
{
  int i, j;
  double A = 1.0, B = 0.0;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
       kdiff[i][j] = ...;
   }
}

void calc_sources(int nx, int ny, double *x, double *y, double *dx, double *dy, double **T, double **b)
{
  int i, j;

  // calculate source (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     b[i][j] = ...;                      // set source function here
     b[i][j] = -b[i][j] * dx[0] * dx[0]; // only works for uniform mesh
   }
}

void set_boundary_conditions(int nx, int ny, double *x, double *y, double *dx, double *dy, double **T, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

  // left boundary -- Homogeneous Dirichlet
  // q(y) * T + p(y) * dTdx = r(y)
  for(j=0; j<ny; j++)
  {
    bcleft[j][0] = 0.0;  // p(y)
    bcleft[j][1] = 1.0;  // q(y)
    bcleft[j][2] = 0.0;  // r(y)
  }

  // right boundary -- Homogeneous Dirichlet
  // q(y) * T + p(y) * dTdx = r(y)
  for(j=0; j<ny; j++)
  {
    bcrght[j][0] = ...;  // p(y)
    bcrght[j][1] = ...;  // q(y)
    bcrght[j][2] = ...;  // r(y)
  }

  // bottom boundary -- Homogeneous  Neumann
  // q(x) * T + p(x) * dTdy = r(x)
  ....;

  // top boundary -- Homogeneous Dirichlet
  // q(x) * T + p(x) * dTdy = r(x)
  ....;
}

void get_coeffs(int nx, int ny, double *x, double *y, double *dx, double *dy, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **kdiff, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

  int i, j;
  double pj, qj, rj;
  double hxhy_sq;

  // calculate diffusivity at [x, y]; may be dependent on T
  calc_diffusivity(nx, ny, x, y, T, kdiff);

  // calculate sources Su, Sp at [x, y]; may be dependent on T
  calc_sources(nx, ny, x, y, dx, dy, T, b);

  // populate values in BC arrays
  set_boundary_conditions(nx, ny, x, y, dx, dy, T, bcleft, bcrght, bctop, bcbot);

  // start populating the coefficients
  hxhy_sq = dx[0]*dx[0]/(dy[0]*dy[0]);

  // ------ Step 1 :: interior points ------
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
   {
     aP[i][j] = ....;
     aE[i][j] = ....;   aW[i][j] = ....;
     aN[i][j] = ....;   aS[i][j] = ....;
   }

  // ------ Step 2 :: left boundary ----------
  i = 0;
   for(j=0; j<ny; j++)
   {
     pj = bcleft[j][0]; qj = bcleft[j][1]; rj = bcleft[j][2];

     aP[i][j] = 2.0 * (1.0 + hxhy_sq) * pj - 2.0*dx[0]*qj;
     aE[i][j] = -2.0*pj;       aW[i][j] = 0.0;
     aN[i][j] = -hxhy_sq * pj;   aS[i][j] = -hxhy_sq * pj;
     b[i][j] = b[i][j] * pj - 2.0 * dx[0] * rj;
   }
  // ------ Step 2 :: left boundary done ---

  // ------ Step 3 :: right boundary ----------
  i = nx-1;
   for(j=0; j<ny; j++)
   {
     pj = bcrght[j][0]; qj = bcrght[j][1]; rj = bcrght[j][2];

     aP[i][j] = ...;
     aW[i][j] = ...;   aE[i][j] = ...;
     aN[i][j] = ...;   aS[i][j] = ...;
     b[i][j] = ....;
   }
  // ------ Step 3 :: right boundary done ---

  // ------ Step 4 :: bottom boundary ----------
  j = 0;
  for(i=1; i<nx-1; i++)
   {
     pj = bcbot[i][0]; qj = bcbot[i][1]; rj = bcbot[i][2];

     aP[i][j] = ...;
     aW[i][j] = ...;   aE[i][j] = ...;
     aN[i][j] = ...;   aS[i][j] = ...;
     b[i][j] = ....;
   }
  // ------ Step 4 :: bottom boundary done ---

  // ------ Step 5 :: top boundary ----------
  j = ny-1;
  for(i=1; i<nx-1; i++)
   {
     pj = bctop[i][0]; qj = bctop[i][1]; rj = bctop[i][2];

     aP[i][j] = ...;
     aW[i][j] = ...;   aE[i][j] = ...;
     aN[i][j] = ...;   aS[i][j] = ...;
     b[i][j] = ....;
   }
  // ------ Step 5 :: top boundary done ---

  // debug
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], b[i][j]);
   }
 
}

double get_max_of_array(int nx, int ny, double **arr)
{
  int i, j;
  double arrmax, val;

  arrmax = arr[0][0];
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = fabs(arr[i][j]);
     if(arrmax < val)
       arrmax = val;
   }
  return arrmax;
}

double get_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{
  double l2err = 0.0, val;
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = arr1[i][j] - arr2[i][j];
     l2err += val*val;
   }
   //printf("l2err = %lf\n", l2err);
  l2err = l2err/((double) (nx*ny));
  l2err = sqrt(l2err);

  return l2err;
}

void solve_gssor(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **Tpad, double **Tpnew, int max_iter, double tol, double relax_T)
{

  int i, j, ip, jp, iter;
  double l2err, arrmax1, arrmax2, err_ref, T_gs, rel_err;

  // copy to padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      Tpad[ip][jp] = T[i][j];
   }

 // now perform iterations
  for(iter=0; iter<max_iter; iter++)
  {
    // update Tpnew
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       T_gs = ( b[i][j] - aE[i][j]*Tpnew[ip+1][jp] - aW[i][j]*Tpnew[ip-1][jp] - 
                          aN[i][j]*Tpnew[ip][jp+1] - aS[i][j]*Tpnew[ip][jp-1] ) / aP[i][j];
       Tpnew[ip][jp] = (1.0-relax_T)*Tpad[ip][jp] + relax_T*T_gs;
       //printf("+++%d %d %e\n", ip, jp, Tpnew[ip][jp]);
     }

    // check for convergence
    l2err = get_l2err_norm(nx+2, ny+2, Tpad, Tpnew);
    arrmax1 = get_max_of_array(nx+2, ny+2, Tpad );
    arrmax2 = get_max_of_array(nx+2, ny+2, Tpnew);
    err_ref = fmax(arrmax1, arrmax2);  err_ref = fmax(err_ref, 1.0e-6);
    rel_err = l2err/err_ref;
    printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
    if(rel_err < tol)
      break;

    // prepare for next iteration
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       Tpad[ip][jp] = Tpnew[ip][jp];
     }

  }

  // copy from padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      T[i][j] = Tpad[ip][jp];
   }

  printf("    >> Fin: %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
}

void solve_jacobi(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **Tpad, double **Tpnew, int max_iter, double tol, double relax_T)
{

  int i, j, ip, jp, iter;
  double l2err, arrmax1, arrmax2, err_ref, T_jac, rel_err;

  // copy to padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      Tpad[ip][jp] = T[i][j];
   }

 // now perform iterations
  for(iter=0; iter<max_iter; iter++)
  {
    // update Tpnew
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       T_jac = ( b[i][j] - aE[i][j]*Tpad[ip+1][jp] - aW[i][j]*Tpad[ip-1][jp] - 
                           aN[i][j]*Tpad[ip][jp+1] - aS[i][j]*Tpad[ip][jp-1] ) / aP[i][j];
       Tpnew[ip][jp] = T_jac;
       //printf("+++%d %d %e\n", ip, jp, Tpnew[ip][jp]);
     }

    // check for convergence
    l2err = get_l2err_norm(nx+2, ny+2, Tpad, Tpnew);
    arrmax1 = get_max_of_array(nx+2, ny+2, Tpad );
    arrmax2 = get_max_of_array(nx+2, ny+2, Tpnew);
    err_ref = fmax(arrmax1, arrmax2);  err_ref = fmax(err_ref, 1.0e-6);
    rel_err = l2err/err_ref;
    printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
    if(rel_err < tol)
      break;

    // prepare for next iteration
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       Tpad[ip][jp] = Tpnew[ip][jp];
     }

  }

  // copy from padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      T[i][j] = Tpad[ip][jp];
   }

  printf("    >> Fin: %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
}

void get_exact_soln(int nx, int ny, double *x, double *y, double **Tex)
{
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     Tex[i][j] = ...;
}

void output_soln(int nx, int ny, int iter, double *x, double *y, double **T, double **Tex)
{
  int i, j;
  FILE* fp;
  char fname[100];

  sprintf(fname, "T_xy_%03d_%03d_%04d.dat", nx, ny, iter);

  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     fprintf(fp, "%lf %lf %lf %lf\n", x[i], y[j], T[i][j], Tex[i][j]);
  fclose(fp);

  printf("Done writing solution for stamp = %d to file %s\n", iter, fname);
}

int main()
{

  int nx, ny;
  double *x, *dx, *y, *dy;
  double ** aP, **aE, **aW, **aN, **aS, **b, **Sp;
  double **T, **Tex, **kdiff, **wrk1, **wrk2;
  double xst, xen, yst, yen;
  double **bcleft, **bcrght, **bctop, **bcbot;
  int i, j, max_iter;
  double relax_T, tol, l2err;
  FILE* fp;  

  // read inputs
  fp = fopen("input.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fclose(fp);

  printf("Inputs are: %d %d %lf %lf %lf %lf\n", nx, ny, xst, xen, yst, yen);

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  x  = (double *)malloc(nx     * sizeof(double));   // grid points
  dx = (double *)malloc((nx-1) * sizeof(double));   // spacing betw grid points

  y  = (double *)malloc(ny     * sizeof(double));   // grid points
  dy = (double *)malloc((ny-1) * sizeof(double));   // spacing betw grid points

  printf("   >> Done allocating 1D arrays -- \n");


  // allocate 2D arrays dynamically
  // -- for T --
  T = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    T[i] = (double *)malloc(ny*sizeof(double));

  // -- for Tex --
  Tex = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    Tex[i] = (double *)malloc(ny*sizeof(double));

  // -- for aP --
  aP = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aP[i] = (double *)malloc(ny*sizeof(double));

  // -- for aE --
  aE = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aE[i] = (double *)malloc(ny*sizeof(double));

  // -- for aW --
  aW = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aW[i] = (double *)malloc(ny*sizeof(double));

  // -- for aN --
  aN = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aN[i] = (double *)malloc(ny*sizeof(double));

  // -- for aS --
  aS = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aS[i] = (double *)malloc(ny*sizeof(double));

  // -- for b --
  b = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    b[i] = (double *)malloc(ny*sizeof(double));

  // left boundary condition
  bcleft  = (double **)malloc( ny * sizeof(double*));
  for(i=0; i<ny; i++)
    bcleft[i] = (double *)malloc(3*sizeof(double));

  // right boundary condition
  bcrght  = (double **)malloc( ny * sizeof(double*));
  for(i=0; i<ny; i++)
    bcrght[i] = (double *)malloc(3*sizeof(double));

  // bottom boundary condition
  bctop   = (double **)malloc( nx * sizeof(double*));
  for(i=0; i<nx; i++)
    bctop[i] = (double *)malloc(3*sizeof(double));

  // top boundary condition
  bcbot   = (double **)malloc( nx * sizeof(double*));
  for(i=0; i<nx; i++)
    bcbot[i] = (double *)malloc(3*sizeof(double));

  // -- for kdiff --
  kdiff = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    kdiff[i] = (double *)malloc(ny*sizeof(double));

  // -- for work arrays wrk1, wrk2 --
  wrk1 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk1[i] = (double *)malloc((ny+2)*sizeof(double));

  wrk2 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk2[i] = (double *)malloc((ny+2)*sizeof(double));


  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");
  
  // initialize the grid
  grid(nx, xst, xen, x, dx);  // -- along x --
  grid(ny, yst, yen, y, dy);  // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  set_initial_guess(nx,ny,x,y,T);  // initial condition
  printf("\n > Done setting up initial guess -- \n");

  // ---
     get_coeffs(nx, ny, x, y, dx, dy,              // grid vars
                aP, aE, aW, aN, aS, b, T, kdiff,   // coefficients
                bcleft, bcrght, bcbot, bctop);     // BC vars
     printf("\n > Done calculating coeffs ----- \n");

     printf("\n > Solving for T ------------- \n");
     max_iter = 100000; tol = 1.0e-10; relax_T = 1.0;
     solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, wrk1, wrk2, max_iter, tol, relax_T);
  // ---
  printf(" > Done solving for T ------------- \n");

  get_exact_soln(nx, ny, x, y, Tex);
  output_soln(nx, ny, 0, x, y, T, Tex);

  l2err = get_l2err_norm(nx, ny, T, Tex);
  printf("%d %d %9.5e", nx, ny, l2err);

  // free memory
   // ----1D arrays ---
   free(y);  free(dy);
   free(x);  free(dx);
   // --- Done 1D arrays ---

   // ----2D arrays ---
   for(i=0; i<nx; i++)     free(T[i]);      free(T);
   for(i=0; i<nx; i++)     free(Tex[i]);    free(Tex);
   for(i=0; i<nx; i++)     free(b[i]);      free(b);
   for(i=0; i<nx; i++)     free(aP[i]);     free(aP);
   for(i=0; i<nx; i++)     free(aE[i]);     free(aE);
   for(i=0; i<nx; i++)     free(aW[i]);     free(aW);
   for(i=0; i<nx; i++)     free(aN[i]);     free(aN);
   for(i=0; i<nx; i++)     free(aS[i]);     free(aS);
   for(i=0; i<nx; i++)     free(kdiff[i]);  free(kdiff);
   for(i=0; i<ny; i++)     free(bcleft[i]); free(bcleft);
   for(i=0; i<ny; i++)     free(bcrght[i]); free(bcrght);
   for(i=0; i<nx; i++)     free(bctop[i]);  free(bctop);
   for(i=0; i<nx; i++)     free(bcbot[i]);  free(bcbot);
   for(i=0; i<nx+2; i++)   free(wrk1[i]);   free(wrk1);
   for(i=0; i<nx+2; i++)   free(wrk2[i]);   free(wrk2);
   // --- Done 2D arrays ---
  printf("\n > Done freeing up memory --------- \n");
  return 0;
}

