#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size;
    int data = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Process 0 initializes the value
        data = 42;
        printf("Process 0 broadcasts data: %d\n", data);
    }

    // All processes receive the broadcasted value
    // Call broadcast function here


    printf("Process %d received data: %d\n", rank, data);

    MPI_Finalize();
    return 0;
}

