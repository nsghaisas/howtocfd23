#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define NX 100         // Number of grid points in x-direction
#define NY 100         // Number of grid points in y-direction
#define L 1.0          // Length of the domain in x-direction
#define W 1.0          // Width of the domain in y-direction
#define TOLERANCE 1e-6 // Tolerance for convergence
#define MAX_ITER 1000  // Maximum number of iterations

double u[NX][NY];        // Solution array
double u_new[NX][NY];    // Updated solution array

int main() {
    int rank, size;
    int i, j, iter;
    double dx, dy;
    double max_diff, global_max_diff;
    double diff;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    dx = L / (NX - 1);
    dy = W / (NY - 1);

    // Initialize the solution
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            u[i][j] = 0.0;
        }
    }

    // Apply boundary conditions
    if (rank == 0) {
        for (j = 0; j < NY; j++) {
            u[0][j] = 1.0;
        }
    }

    if (rank == size - 1) {
        for (j = 0; j < NY; j++) {
            u[NX - 1][j] = 1.0;
        }
    }

    // Perform iterations until convergence or maximum iterations reached
    for (iter = 0; iter < MAX_ITER; iter++) {
        max_diff = 0.0;

        // Update the interior points
        for (i = 1; i < NX - 1; i++) {
            for (j = 1; j < NY - 1; j++) {
                u_new[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                diff = fabs(u_new[i][j] - u[i][j]);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
        }

        // Exchange boundary values between neighboring processes : TODO



        // Compute the maximum difference across all processes : TODO



        // Copy the updated solution to the solution array
        for (i = 1; i < NX - 1; i++) {
            for (j = 0; j < NY; j++) {
                u[i][j] = u_new[i][j];
            }
        }

        // Check for convergence
        if (global_max_diff < TOLERANCE) {
            break;
        }
    }

    // Print the solution
    if (rank == 0) {
        for (i = 0; i < NX; i++) {
            for (j = 0; j < NY; j++) {
                printf("%.4f ", u[i][j]);
            }
            printf("\n");
        }
    }

    MPI_Finalize();
    return 0;
}

