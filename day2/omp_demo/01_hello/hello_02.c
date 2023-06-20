/* use function call
 * int tid */

#include<stdio.h>
#include<omp.h>

int main()
{
  int p, i;
  p = 24;
  omp_set_num_threads(p);

  #pragma omp parallel
  {
    int tid; 
    tid = omp_get_thread_num();
    printf("Hi from thread number %d\n", tid);
  }

  return 0;
}
