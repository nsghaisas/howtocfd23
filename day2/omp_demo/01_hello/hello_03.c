/* use function call to set number of threads
 * int tid
 * function call from parallel region */

#include<stdio.h>
#include<omp.h>

void work(int tid)
{
  printf("Hi from a function called from inside a parallel region %d\n", tid);
}

int main()
{
  int p, i;
  p = 4;
  omp_set_num_threads(p);

  #pragma omp parallel
  {
    int tid; 
    tid = omp_get_thread_num();
    //printf("Hi from thread number %d\n", tid);
    work(tid);
  }

  return 0;
}
