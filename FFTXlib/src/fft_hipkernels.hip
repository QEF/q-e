#include <hip/hip_runtime.h>
#include <stdio.h>

__global__ void scalar_multiply_kernel(double *dev_ptr, double val, int s) 
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
        int stride = gridDim.x * blockDim.x;
	for (; i < s; i += stride) {
		dev_ptr[i] *= val;
	}
}

extern "C" { 
  void c_scalar_multiply_(int s, double val, double *dev_ptr, hipStream_t stream)
  {
   
    const int blocks = 32;
    const int threads = 256;
    hipError_t err;

    hipLaunchKernelGGL(scalar_multiply_kernel, dim3(blocks), dim3(threads), 0, stream, dev_ptr, val, s );
  }
}

__global__ void loop2d_scatter_kernel_in(double *dev_ptr2, double *dev_ptr, int *pmap, unsigned nppx, unsigned nnp, unsigned of2, int of, unsigned npp, unsigned nswip)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x+1;
        int j = blockIdx.y * blockDim.y + threadIdx.y+1;

        if (j<=npp && i<=nswip) {
           dev_ptr[pmap[i-1]*2+of+(j-1)*nnp*2 -2] = dev_ptr2[j*2+of2+(i-1)*nppx*2 -2];
           dev_ptr[pmap[i-1]*2+of+(j-1)*nnp*2 -1] = dev_ptr2[j*2+of2+(i-1)*nppx*2 -1];
        }
}

__global__ void loop2d_scatter_kernel_fw(double *dev_ptr, double *dev_ptr2, int *pmap, unsigned nppx, unsigned nnp, unsigned of2, int of, unsigned npp, unsigned nswip)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x+1;
        int j = blockIdx.y * blockDim.y + threadIdx.y+1;

        if (j<=npp && i<=nswip) {
           dev_ptr2[j*2+of2+(i-1)*nppx*2 -2] = dev_ptr[pmap[i-1]*2+of+(j-1)*nnp*2 -2];
           dev_ptr2[j*2+of2+(i-1)*nppx*2 -1] = dev_ptr[pmap[i-1]*2+of+(j-1)*nnp*2 -1];
        }
}


extern "C" {
  void loop2d_scatter_hip_(int drz, double *dev_ptr2, double *dev_ptr, int *pmap, unsigned nppx, unsigned nnp, unsigned of2, int of, unsigned npp, unsigned nswip, hipStream_t stream ) {

//    const int blocks = 32;
//   const int threads = 256;
 
   const int blocks_i = (nswip+31)/32;
   const int blocks_j = (npp+31)/32;
   const int threads_i = 32;
   const int threads_j = 32;

   hipError_t err;

   if (drz==-1) {
      hipLaunchKernelGGL(loop2d_scatter_kernel_in, dim3(blocks_i,blocks_j), dim3(threads_i,threads_j), 0, stream, dev_ptr2, dev_ptr, pmap, nppx, nnp, of2, of, npp, nswip );
   }
   if (drz== 1) {
      hipLaunchKernelGGL(loop2d_scatter_kernel_fw, dim3(blocks_i,blocks_j), dim3(threads_i,threads_j), 0, stream, dev_ptr2, dev_ptr, pmap, nppx, nnp, of2, of, npp, nswip );
   }
//   err = hipDeviceSynchronize();
  }
}

