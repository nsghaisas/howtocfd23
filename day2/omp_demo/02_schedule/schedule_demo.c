/* deomstrate how work is split by schedule clause and options
 * default
 * static  (4|8|16|32)
 * dynamic (1|4|8)
 * guided  (1|2|4) */

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

int main()
{
  int p = 4, N = 64, *thread_assigned;
  int i, ip;

  // allocate and initialize to -1
  thread_assigned = (int*) malloc(N*sizeof(int));
  for(i=0; i<N; i++)
    thread_assigned[i] = -1;


  // set number of threads
  omp_set_num_threads(p);

  #pragma omp parallel
  {
    #pragma omp for schedule (dynamic, 4)
    for(i=0; i<N; i++)
    {
      int tid = omp_get_thread_num();
      thread_assigned[i] = tid;
    }
  }

  // visualize thread_assigned
  for(ip=0; ip<p; ip++)
  {
    printf("\n");
    for(i=0; i<N; i++)
    {
      if(thread_assigned[i]==ip)
        printf("X");
      else
        printf(" ");
    }
  }
  printf("\n");

  return 0;
}
