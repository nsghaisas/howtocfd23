#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size;
    int my_value = 0;
    int sum = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    my_value = rank + 1;  // Assign a unique value to each process

    // Perform the reduction operation across all processes
    // Call AllReduce function here: TODO

  
    printf("Process %d: my_value = %d, sum = %d\n", rank, my_value, sum);

    MPI_Finalize();
    return 0;
}

