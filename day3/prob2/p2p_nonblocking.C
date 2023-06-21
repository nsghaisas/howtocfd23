#include <stdio.h>
#include <mpi.h>

int main() {
    int rank, size;
    int data_send = 42;
    int data_recv;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        fprintf(stderr, "This program requires at least 2 processes.\n");
        MPI_Finalize();
        return 1;
    }

    MPI_Request request_send, request_recv;
    MPI_Status status;

    if (rank == 0) {
        // Process 0 sends a message to process 1
	// Call non-blocking send function here : TODO


        printf("Process 0 initiated non-blocking send of data: %d\n", data_send);
    } else if (rank == 1) {
        // Process 1 receives the message from process 0
	// Call non-blocking recv function here : TODO


        printf("Process 1 initiated non-blocking receive\n");

        // Wait for the receive operation to complete
        // Call wait function here : TODO

        printf("Process 1 received data: %d\n", data_recv);
    }

    MPI_Finalize();
    return 0;
}

