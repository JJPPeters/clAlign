// This class is basically just a container to easily pass around all the neede kernels etc

// TODO: more protection on accessing kernels etc if no device has been set (i.e. getters with error checking.)

#pragma once

#include "stdafx.h"
#include "clWrapper/clWrapper.h"
#include "DMWrapper/Dmout.h"
#include "boost/thread/mutex.hpp"

#include "Kernels.h"

typedef boost::shared_ptr<clContext> shared_context;
typedef boost::shared_ptr<clFourier> shared_fft;
typedef boost::shared_ptr<clKernel> shared_kernel;

// Not really OpenCL but an easy place to have this for now
typedef boost::shared_ptr<boost::mutex> shared_mutex;

class clArgStore
{
public:
	shared_context Context;

	shared_fft FFT;

	shared_kernel kMultiCorrelation;
	shared_kernel kFFTShift;
	shared_kernel kBilinearInterpolate;
	shared_kernel kExponentialPass;

	bool haveDevice;

	clArgStore()  : haveDevice(false) {}

	clArgStore(clDevice dev) : haveDevice(false) { SetContext(dev); }

	void SetContext(clDevice dev)
	{
		haveDevice = false;
		Context = boost::make_shared<clContext>(OpenCL::MakeContext(dev, Queue::InOrder));

		kMultiCorrelation.reset(new clKernel(*Context, sMultiCorrelation, 6, "clMultiCorrelation"));
		kExponentialPass.reset(new clKernel(*Context, sExponentialPass, 4, "clExponentialPass"));
		kFFTShift.reset(new clKernel(*Context, sfftShift, 4, "clfftShift"));
		kBilinearInterpolate.reset(new clKernel(*Context, sBilinearInterpolate, 14, "clBilinearInterpolate"));
		
		if (CheckStatus("ERROR: Unable to set up OpenCL with code "))
			haveDevice = true;
	};

	bool CheckStatus(std::string message)
	{
		bool ok = Context->GetStatus() == CL_SUCCESS;
		if (!ok)
			DMresult << "Error: " << message << ". Error code " << Context->GetStatus() << DMendl;
		return ok;
	};
};