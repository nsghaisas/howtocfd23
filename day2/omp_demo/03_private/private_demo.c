/* deomstrate private, firstprivate, lastprivate, reduction*/

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

int main()
{
  int p = 4; //, N = 64, *thread_assigned;
  int i, j; //, ip;


  // set number of threads
  omp_set_num_threads(p);

  i = 10;

  #pragma omp parallel firstprivate(i) 
  {
    int tid = omp_get_thread_num();
    printf("Before update, tid = %d, i = %d\n", tid, i);
    i += (tid+1)*1000;
    printf("After  update, tid = %d, i = %d\n", tid, i);
  }
  printf("After parallel region, i = %d\n", i);

    // Also try:
    //#pragma omp parallel firstprivate(i)
    //#pragma omp parallel private(i)
    //#pragma omp parallel



  //#pragma omp parallel
  //{ 
  //  int tid = omp_get_thread_num();
  //  #pragma omp sections firstprivate(i) lastprivate(i)
  //  {
  //    #pragma omp section
  //    {
  //      printf("Section 1 before update, tid = %d, i = %d\n", tid, i);
  //      i += (tid+1)*1000;
  //      printf("Section 1 after  update, tid = %d, i = %d\n", tid, i);
  //      
  //    }
  //    #pragma omp section
  //    {
  //      printf("Section 2 before update, tid = %d, i = %d\n", tid, i);
  //      i += (tid+1)*1000;
  //      printf("Section 2 after  update, tid = %d, i = %d\n", tid, i);
  //    }
  //    #pragma omp section
  //    {
  //      printf("Section 3 before update, tid = %d, i = %d\n", tid, i);
  //      i += (tid+1)*1000;
  //      printf("Section 3 after  update, tid = %d, i = %d\n", tid, i);
  //    }
  //  }
  //}
  //printf("After parallel sections, i = %d\n", i);


    // Also try:
    //#pragma omp sections private(i)
    //#pragma omp sections firstprivate(i)
    //#pragma omp sections firstprivate(i) lastprivate(i)

  return 0;
}
