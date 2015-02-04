#include "stdafx.h"
#include "clArgStore.h"

shared_context clArgStore::Context;

shared_fft clArgStore::FFT;

shared_kernel clArgStore::kMultiCorrelation;
shared_kernel clArgStore::kFFTShift;
shared_kernel clArgStore::kBilinearInterpolate;
shared_kernel clArgStore::kExponentialPass;

std::vector<shared_complex> clArgStore::ComplexBuffers(4);

bool clArgStore::haveDevice;