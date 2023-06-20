#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<omp.h>

int main()
{

  double *a, *local_sum, p = 8;
  int i, N = 1000;
  double sum, sum_reduction;

  omp_set_num_threads(p);

  a = (double*) malloc(N*sizeof(double));
  if(a==NULL) { printf("Could not allocate a. Please check.\n"); exit(0); }

  local_sum = (double*) malloc(p*sizeof(double));
  if(local_sum==NULL) { printf("Could not allocate local_sum. Please check.\n"); exit(0); }


  // initialize a
  for (i=0; i<N; i++)
  {
    a[i] = (double)i;
  }

  //-------------------------------------------------------
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    #pragma omp for
    for(i=0; i<N; i++)
    {
      local_sum[tid] += a[i];
    }
  }

  // compile all local_sums into sum and output
  sum = 0.0;
  for(i=0; i<p; i++)
    sum += local_sum[i];
  printf("Sum of all numbers = %f\n", sum);
  //-------------------------------------------------------


  // Alternative: use reduction
  sum_reduction = 0.0;
  #pragma omp parallel
    #pragma omp for reduction(+:sum_reduction)
    for (i=0; i<N; i++)
      sum_reduction += a[i];
  printf("Sum with reduction = %f\n", sum_reduction);

  if(a!=NULL) free(a);
  if(local_sum!=NULL) free(local_sum);

  return 0;
}
