module load gnu8
gcc sin_fft3.c -I/hpc/apps/PadeOpsDepGNU/fftw-3.3.5/include -L/hpc/apps/PadeOpsDepGNU/fftw-3.3.5/lib -lfftw3 -lm
#gfortran sin_fft3.f90 -I/hpc/apps/PadeOpsDepGNU/fftw-3.3.5/include -L/hpc/apps/PadeOpsDepGNU/fftw-3.3.5/lib -lfftw3 -lm
