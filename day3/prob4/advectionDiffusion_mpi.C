#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846264338327
#define NX 100       // Number of grid points
#define NT 1000      // Number of time steps
#define L 2 * PI        // Length of the domain
#define T_FINAL 1.0  // Final simulation time
#define D 0.1        // Diffusion coefficient
#define U 1.0        // Advection velocity
#define DT (T_FINAL / NT)         // Time step size
#define DX (L / (NX))             // Grid spacing
#define ALPHA (D * DT / (DX * DX)) // Stability parameter

int main(int argc, char** argv) {
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int chunk_size = NX / num_procs;
    int start_index = rank * chunk_size;
    int end_index = start_index + chunk_size;

    if (rank == num_procs - 1) {
        // Adjust the end_index for the last process to account for any remaining grid points
        end_index = NX;
    }

    double u[chunk_size + 2];         // Array to store the solution for each process (with ghost cells)
    double u_new[chunk_size + 2];     // Array to store the updated solution for each process (with ghost cells)
    double x;                        // Current position
    double t;                        // Current time
    int i, n;                        // Loop variables

    // Initialize the solution for each process (including ghost cells)
    for (i = 1; i <= chunk_size; i++) {
        x = (i + start_index - 1) * DX;
        u[i] = sin(x);
    }

    // Apply periodic boundary conditions
    u[0] = u[chunk_size];
    u[chunk_size + 1] = u[1];

    // Perform time integration
    for (n = 0; n < NT; n++) {
        t = n * DT;

        // Update the solution using Euler integration and central difference
        for (i = 1; i <= chunk_size; i++) {
            u_new[i] = u[i] + ALPHA * (u[i + 1] - 2 * u[i] + u[i - 1]) +
                       U * DT * (u[i + 1] - u[i - 1]) / (2 * DX);
        }

        // Exchange ghost cell data between neighboring processes : TODO



        // Update the solution array (excluding ghost cells)
        for (i = 1; i <= chunk_size; i++) {
            u[i] = u_new[i];
        }
    }

    // Gather the final solution from all processes to process 0
    double* final_u = NULL;
    if (rank == 0) {
        final_u = (double*)malloc(NX * sizeof(double));
    }

    MPI_Gather(&u[1], chunk_size, MPI_DOUBLE, final_u, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print the final solution (only done by process 0)
    if (rank == 0) {
        for (i = 0; i < NX; i++) {
            x = i * DX;
            printf("x = %.2f\tu = %.4f\n", x, final_u[i]);
        }

        free(final_u);
    }

    MPI_Finalize();
    return 0;
}

