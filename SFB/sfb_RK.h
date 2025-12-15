#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cublas_v2.h>

#include "types.h"

#include "GPUkernels_RK.cu"
#include "magma_v2.h"
#include "magma_lapack.h"

#include "input.h"
#include "EulerQuaternion.h"

