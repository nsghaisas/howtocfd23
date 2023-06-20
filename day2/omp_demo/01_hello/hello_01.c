/* Check the number of cores on the system
 * sse, sse2
 * set OMP_NUM_THREADS in terminal/ bashrc 
 * use function call
 * clause to parallel construct
 * thread id before and after parallel  */

#include<stdio.h>
#include<omp.h>

int main()
{
  int p, i;
  p = 4;
  //omp_set_num_threads(p);

  #pragma omp parallel
  {
    //printf("Hi from thread number %d\n",omp_get_thread_num());
    printf("Hi from a parallel region\n");
  }

  return 0;
}
