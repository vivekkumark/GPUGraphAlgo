#ifndef LSG_COMMON
#define LSG_COMMON

#include <stdio.h>
#include <cuda.h>
#include <time.h>
#include <fstream>
#include <string>
#include <iostream>
#include <limits>
#include <string.h>

#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <inttypes.h>

//#define DATA_DRIVEN //Partially implemented

//First stage : without using 3-phase technique
//#define ENABLE_ATOMIC_MARKING

#define ENABLE_SECONDARY_MARKING 1

//extra marking for increasing the sunset cavity(But observed very poor performance)
//#define EXTRA_MARKIN 0

//#define ENABLE_POST_WORK
//#define ENABLE_DYNAMICP

#define CAVITY_SIZE 150
#define BLOCKSIZE 512
#define _BS 768

typedef unsigned foru;
using namespace std;

static unsigned CudaTest(char *msg)
{
  cudaError_t e;
  cudaDeviceSynchronize();
  if (cudaSuccess != (e = cudaGetLastError())) {
    fprintf(stderr, "%s: %d\n", msg, e);
    fprintf(stderr, "%s\n", cudaGetErrorString(e));
    //exit(-1);
    return 1;
  }
  return 0;
}

#endif
