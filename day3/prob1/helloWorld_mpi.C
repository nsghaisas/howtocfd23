#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment : TODO

    // Get the number of processes : TODO
    int npes;



    // Get the rank of the process : TODO
    int myid;

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello World! from processor %s, rank %d out of %d processors\n",
           processor_name, myid, npes);

    // Finalize the MPI environment : TODO 
}

