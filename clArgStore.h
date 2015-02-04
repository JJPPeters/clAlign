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
typedef boost::shared_ptr<clMemory<std::complex<float>, Auto>> shared_complex;


// Not really OpenCL but an easy place to have this for now
typedef boost::shared_ptr<boost::mutex> shared_mutex;

static class clArgStore
{
public:
	static shared_context Context;

	static shared_fft FFT;

	static shared_kernel kMultiCorrelation;
	static shared_kernel kFFTShift;
	static shared_kernel kBilinearInterpolate;
	static shared_kernel kExponentialPass;

	static std::vector<shared_complex> ComplexBuffers;

	static bool haveDevice;

	clArgStore() { haveDevice = false; }

	//clArgStore(clDevice dev) : haveDevice(false) { SetContext(dev); }

	//~clArgStore(){ DMresult << "Destructing argstore" << DMendl; }

	static void SetContext(clDevice dev)
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

	static bool CheckStatus(std::string message)
	{
		bool ok = Context->GetStatus() == CL_SUCCESS;
		if (!ok)
			DMresult << "Error: " << message << ". Error code " << Context->GetStatus() << DMendl;
		return ok;
	};

	static bool CreateComplex(int width, int height)
	{
		for (int i = 0; i < ComplexBuffers.size(); i++)
		{
			DMresult << "Before " << i << " ref count = " << ComplexBuffers[i].use_count() << DMendl;
			ComplexBuffers[i].reset();
			DMresult << "Before " << i << " ref count = " << ComplexBuffers[i].use_count() << DMendl;
			ComplexBuffers[i] = (*Context).CreateBuffer<std::complex<float>, Auto>(width * height);
			DMresult << "Middle " << i << DMendl;
			if (!CheckStatus("Creating buffer")) return false;
			DMresult << "After " << i << DMendl;
		}
		return true;
	}
};