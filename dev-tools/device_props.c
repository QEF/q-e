/*
 * device_props
 *
 * Copyright (c) 2017 Pietro Bonfa', Adapted from Simatra Modelling Technologies
 *
 * Returns a sorted list of available CUDA devices with properties colon delimited.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

#include <cuda_runtime_api.h>

#if defined( __APPLE__ ) || defined(  __MACH__ )
#define CUDART_LIBRARY_NAME "libcudart.dylib"
#elif defined( unix ) || defined( __unix ) || defined( __unix__ ) || defined( __linux__ ) || defined( __FreeBSD__ )
#define CUDART_LIBRARY_NAME "libcudart.so"
#elif defined(_WIN32)
#error Windows not implemented.
#else
#error Must define CUDART_LIBRARY_NAME (e.g. libcudart.so or libcudart.dylib)
#endif

enum status {
  DeviceProps_Success,
  DeviceProps_EmulationOnly,
  DeviceProps_NoDevices,
  DeviceProps_UnknownError,
  DeviceProps_NoCudaRuntime,
  DeviceProps_NoRuntimeVersion,
  DeviceProps_NoDriverVersion
};


typedef struct{
  struct cudaDeviceProp props;
  int gflops;
  int unsorted;
}simCudaDevice;

// Function pointer types to dynamically loaded functions from libcudart
typedef cudaError_t (*cudaGetDeviceCount_f)(int *);
typedef cudaError_t (*cudaGetDeviceProperties_f)(struct cudaDeviceProp*, int);


#define MAX_DEVICES 20
#define BUFFER_LENGTH 1024*1024
static char props_buffer[MAX_DEVICES * BUFFER_LENGTH];
static char error_message[BUFFER_LENGTH];

char *devicePropsError (void) {
  return error_message;
}

int getSystemProps (char **hostProps) {
  int cudaRuntimeVersion;
  int cudaDriverVersion;
  // Retrive system properties
  if (cudaSuccess != cudaRuntimeGetVersion( &cudaRuntimeVersion )) {
    snprintf(error_message, BUFFER_LENGTH,
	     "Could not retrive Runtime Version.\n"
	     "\tSomething is wrong with your installation.");
    error_message[BUFFER_LENGTH - 1] = '\0';
    return DeviceProps_NoRuntimeVersion;
  }
  if (cudaSuccess != cudaDriverGetVersion( &cudaDriverVersion )) {
    snprintf(error_message, BUFFER_LENGTH,
	     "Could not retrive Runtime Version.\n"
	     "\tSomething is wrong with your installation.");
    error_message[BUFFER_LENGTH - 1] = '\0';
    return DeviceProps_NoDriverVersion;
  }
  

  // Write global info
  char *write = props_buffer;
  write += 1 + snprintf(write, BUFFER_LENGTH, 
                          " driverVersion: %d\n runtimeVersion: %d\n", 
                          cudaDriverVersion, cudaRuntimeVersion);
  *hostProps = props_buffer;
  error_message[0] = '\0';
  return DeviceProps_Success;
}

int getDeviceProps (int *deviceCount, char **deviceProps) {
  // Cuda Runtime interface
  void *cudaRT = NULL;
  cudaGetDeviceCount_f cudaGetDeviceCount = NULL;
  cudaGetDeviceProperties_f cudaGetDeviceProperties = NULL;

  cudaError_t cuErr;
  int ndevices; // Number of devices reported by Cuda runtime
  int undevices = 0; // Number of devices that are unusable by simEngine
  unsigned int deviceid;
  unsigned int sort;
  simCudaDevice *devices;

  cudaRT = dlopen(CUDART_LIBRARY_NAME, RTLD_NOW);
  if(!cudaRT) {
    char full_library_name[PATH_MAX];
    sprintf(full_library_name, "/usr/local/cuda/lib64/%s", CUDART_LIBRARY_NAME);
    cudaRT = dlopen(full_library_name, RTLD_NOW);
    if(!cudaRT) {
      sprintf(full_library_name, "/usr/local/cuda/lib/%s", CUDART_LIBRARY_NAME);
      cudaRT = dlopen(full_library_name, RTLD_NOW);
      if(!cudaRT) {
	snprintf(error_message, BUFFER_LENGTH,
		 "Failed to load CUDA runtime environment from %s.\n"
		 "\tIs the CUDA runtime environment installed in the default location\n"
		 "\tOR is LD_LIBRARY_PATH environment variable set to include CUDA libraries?",
		 CUDART_LIBRARY_NAME);
	error_message[BUFFER_LENGTH - 1] = '\0';
	return DeviceProps_NoCudaRuntime;
      }
    }
  }

  cudaGetDeviceCount = (cudaGetDeviceCount_f)dlsym(cudaRT, "cudaGetDeviceCount");
  cudaGetDeviceProperties = (cudaGetDeviceProperties_f)dlsym(cudaRT, "cudaGetDeviceProperties");

  if(!cudaGetDeviceCount || !cudaGetDeviceProperties) {
    snprintf(error_message, BUFFER_LENGTH, 
	     "Failed to load CUDA functions from %s.\n"
	     "\tThe CUDA library found is incompatible with simEngine.",
	     CUDART_LIBRARY_NAME);
    error_message[BUFFER_LENGTH - 1] = '\0';
    return DeviceProps_NoCudaRuntime;
  }

  if (cudaSuccess != cudaGetDeviceCount(&ndevices)) {
    snprintf(error_message, BUFFER_LENGTH,
	     "Error obtaining device count.\n"
	     "\tIs there a CUDA capable GPU available on this computer?");
    error_message[BUFFER_LENGTH - 1] = '\0';
    return DeviceProps_UnknownError;
  }

  if (0 == ndevices) {
    snprintf(error_message, BUFFER_LENGTH,
	     "No suitable devices found.\n"
	     "\tIs your CUDA driver installed, and have you rebooted since installation?");
    error_message[BUFFER_LENGTH - 1] = '\0';
    return DeviceProps_NoDevices;
  }

  devices = (simCudaDevice *)malloc(sizeof(simCudaDevice) * ndevices);

  // Retrieve the properties for all Cuda devices
  for (deviceid = 0; deviceid < ndevices; ++deviceid) {
    if (cudaSuccess != cudaGetDeviceProperties(&devices[deviceid-undevices].props, deviceid)) {
      snprintf(error_message, BUFFER_LENGTH,
	       "Error obtaining properties for device %d.\n"
	       "\tThe CUDA library found is incompatible with simEngine.", 
	       deviceid);
      error_message[BUFFER_LENGTH - 1] = '\0';
      free(devices);
      return DeviceProps_UnknownError;
    }
    // Filter out emulation devices
    if(9999 == devices[deviceid-undevices].props.major) {
      undevices += 1;
    }
    // Track GFLOPs of real devices
    else {
      devices[deviceid-undevices].gflops = devices[deviceid-undevices].props.multiProcessorCount * devices[deviceid-undevices].props.clockRate;
      devices[deviceid-undevices].unsorted = 1;
    }
  }

  // Subtract emulation devices from device count
  *deviceCount = ndevices - undevices;
  if (0 == *deviceCount) {
    snprintf(error_message, BUFFER_LENGTH,
	     "Only emulation device found.\n"
	     "\tDo you have a CUDA device?\n"
	     "\tIs the CUDA driver installed?\n"
	     "\tHave you rebooted after installing the driver?\n"
	     "\tDo you have device permissions set to allow CUDA computation?");
    error_message[BUFFER_LENGTH - 1] = '\0';
    free(devices);
    return DeviceProps_EmulationOnly;
  }

  // Sort the useable devices by max GFLOPs
  char *write = props_buffer;
  for(sort = 0; sort<(ndevices - undevices) && sort<MAX_DEVICES; ++sort){
    int max_gflops = 0;
    int max_gflops_dev = 0;
    int written = 0;
    for(deviceid = 0; deviceid<(ndevices - undevices); ++deviceid){
      if(devices[deviceid].unsorted && devices[deviceid].gflops > max_gflops){
        max_gflops = devices[deviceid].gflops;
        max_gflops_dev = deviceid;
      }
    }
    // Print one device per line with properties colon separated
    written = snprintf(write, BUFFER_LENGTH,
           // One line output: "%d:%s:%zd:%zd:%d:%d:%zd:%d:%d,%d,%d:%d,%d,%d:%zd:%d:%d:%d:%zd:%d:%d:%d:%d:%d:%d",
           "-\n"
           " devId:  %d\n"
           " name: \"%s\"\n"
           " totalGlobalMem: %zd\n"
           " sharedMemPerBlock: %zd\n"
           " regsPerBlock: %d\n"
           " warpSize: %d\n"
           " memPitch: %zd\n"
           " maxThreadsPerBlock: %d\n"
           " maxThreadsDim[0]: %d\n"
           " maxThreadsDim[1]: %d\n"
           " maxThreadsDim[2]: %d\n"
           " maxGridSize[0]: %d\n"
           " maxGridSize[1]: %d\n"
           " maxGridSize[2]: %d\n"
           " totalConstMem: %zd\n"
           " major: %d\n"
           " minor: %d\n"
           " clockRate: %d\n"
           " textureAlignment: %zd\n"
           " deviceOverlap: %d\n"
           " multiProcessorCount: %d\n"
           " kernelExecTimeoutEnabled: %d\n"
           " integrated: %d\n"
           " canMapHostMemory: %d\n"
           " computeMode: %d\n",
		       max_gflops_dev,
		       devices[max_gflops_dev].props.name,
		       // Switch to kB to not overflow an int
		       devices[max_gflops_dev].props.totalGlobalMem>>10,
		       devices[max_gflops_dev].props.sharedMemPerBlock,
		       devices[max_gflops_dev].props.regsPerBlock,
		       devices[max_gflops_dev].props.warpSize,
		       devices[max_gflops_dev].props.memPitch,
		       devices[max_gflops_dev].props.maxThreadsPerBlock,
		       devices[max_gflops_dev].props.maxThreadsDim[0],
		       devices[max_gflops_dev].props.maxThreadsDim[1],
		       devices[max_gflops_dev].props.maxThreadsDim[2],
		       devices[max_gflops_dev].props.maxGridSize[0],
		       devices[max_gflops_dev].props.maxGridSize[1],
		       devices[max_gflops_dev].props.maxGridSize[2],
		       devices[max_gflops_dev].props.totalConstMem,
		       devices[max_gflops_dev].props.major,
		       devices[max_gflops_dev].props.minor,
		       devices[max_gflops_dev].props.clockRate,
		       devices[max_gflops_dev].props.textureAlignment,
		       devices[max_gflops_dev].props.deviceOverlap,
		       devices[max_gflops_dev].props.multiProcessorCount,
		       devices[max_gflops_dev].props.kernelExecTimeoutEnabled,
		       devices[max_gflops_dev].props.integrated,
		       devices[max_gflops_dev].props.canMapHostMemory,
		       devices[max_gflops_dev].props.computeMode
		       );
    write += 1 + written;
    devices[max_gflops_dev].unsorted = 0;
  }

  *deviceProps = props_buffer;

  free(devices);
  error_message[0] = '\0';
  return DeviceProps_Success;
}


int main(int argc, char **argv){
  int ndevices, sort;
  char *props;
  
  int status = getSystemProps(&props);
  if (DeviceProps_Success != status) {
    fprintf(stderr, "%s", devicePropsError());
    fprintf(stderr, "\n");
    return status;
  }
  // Print global properties
  fprintf(stdout, "system:\n");
  fprintf(stdout, props);
  
  status = getDeviceProps(&ndevices, &props);
  if (DeviceProps_Success != status) {
    fprintf(stderr, "%s", devicePropsError());
    fprintf(stderr, "\n");
    return status;
  }
  // Print device properties
  fprintf(stdout, "devices:\n");
  for (sort = 0; sort < ndevices; sort++) {
    int written = fprintf(stdout, props);
    fprintf(stdout, "\n");
    props += 1 + written;
  }

  return status;
}
