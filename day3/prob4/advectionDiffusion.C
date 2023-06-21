#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846264338327
#define NX 100       // Number of grid points
#define NT 1000      // Number of time steps
#define L 2*PI        // Length of the domain
#define T_FINAL 1.0  // Final simulation time
#define D 0.1        // Diffusion coefficient
#define U 1.0        // Advection velocity
#define DT (T_FINAL / NT)         // Time step size
#define DX (L / (NX))             // Grid spacing
#define ALPHA (D * DT / (DX * DX)) // Stability parameter

int main() {
    double u[NX];         // Array to store the solution
    double u_new[NX];     // Array to store the updated solution
    double x;               // Current position
    double t;               // Current time
    int i, n;               // Loop variables

    // Initialize the solution
    for (i = 1; i < NX - 1; i++) {
        x = (i - 1) * DX;
        u[i] = sin(x);
    }

    // Apply boundary conditions
    u[0] = u[NX - 2];
    u[NX - 1] = u[1];

    // Perform time integration
    for (n = 0; n < NT; n++) {
        t = n * DT;

        // Update the solution using Euler integration and central difference
        for (i = 1; i < NX - 1; i++) {
            u_new[i] = u[i] + ALPHA * (u[i + 1] - 2 * u[i] + u[i - 1]) +
                         U * DT * (u[i + 1] - u[i - 1]) / (2 * DX);
        }

        // Apply boundary conditions
        u_new[0] = u_new[NX - 2];
        u_new[NX - 1] = u_new[1];

        // Update the solution array
        for (i = 0; i < NX; i++) {
            u[i] = u_new[i];
        }
    }

    // Print the final solution
    for (i = 0; i < NX; i++) {
        x = i * DX;
        printf("x = %.2f\tphi = %.4f\n", x,u[i]);
    }

    return 0;
}

