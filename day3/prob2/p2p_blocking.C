#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size;
    int data = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        fprintf(stderr, "This program requires at least 2 processes.\n");
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        // Process 0 sends a message to process 1
        data = 42;

        // Call send function here : TODO

        printf("Process 0 sent data: %d\n", data);
    } else if (rank == 1) {
        // Process 1 receives the message from process 0
        // Call recv function here : TODO

        printf("Process 1 received data: %d\n", data);
    }

    MPI_Finalize();
    return 0;
}

