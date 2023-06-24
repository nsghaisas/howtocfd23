#nvfortran -ta=tesla:cc80 -Minfo=accel -cpp -O3 diffusion_gpu_new.cuf nvtx.F90 -o diff_gpu.exe  -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5//cuda/11.7/lib64/  -lnvToolsExt
##nvfortran -ta=tesla:cc70 -Minfo=accel -cpp -O3 diffusion_openacc_singproc.f90 nvtx.F90 -o diff_gpu_acc.exe  -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5//cuda/11.7/lib64/  -lnvToolsExt
nvfortran -ta=tesla:cc70 -Minfo=accel -cpp -O3 diffusion_openacc_singproc.f90 -o diff_gpu_acc.exe  -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5//cuda/11.7/lib64/  
#nvfortran -O3 diffusion_openacc_singproc.f90 nvtx.F90 -o diff_sp.exe  -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5//cuda/11.7/lib64/  -lnvToolsExt
#mpif90 -O3 diffusion.f90 nvtx.F90 -o diff_mpi.exe  -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5//cuda/11.7/lib64/  -lnvToolsExt
