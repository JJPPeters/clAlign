#pragma once

#include "WorkerClass.h"
#include "clWrapper/clWrapper.h"
#include "clArgStore.h"

#include "DMWrapper/DMImage.h"
#include "DMWrapper/DMColmap.h"
#include "DMWrapper/DMROI.h"

#include <boost/shared_ptr.hpp>
#include "boost/thread/mutex.hpp"
#include "boost/thread.hpp"

class Alignment : public WorkerClass, public DMListenable
{
private:
	// Main body of the alignment code. Runs in new thread to avoid locking anything else.
	void DoWork();

	// Performs an alignment by solving an overdetermined set of equations relating all shifts
	void OverDeterminedAlign();

	// Performs Cross correlation and phase correlation using OpenCL
	std::vector<std::complex<float>> CrossCorrelation(std::vector<std::complex<float>> image1, std::vector<std::complex<float>>image2);

	// Finds the pixel maxima of a vector
	coord<int> FindMaxima(std::vector<std::complex<float>> data);

	// Fit parabola to find refined maxima
	coord<float>FindVertexParabola(std::vector<float> data);

	void AlignImage(std::vector<float> shiftx, std::vector<float> shifty);

	boost::shared_ptr<boost::mutex> _mtx;
	boost::mutex WorkLock;

	DMImage Image;

	int width, height, depth;

	boost::shared_ptr<clArgStore> clArguments;
	std::vector<boost::shared_ptr<clMemory<std::complex<float>, Auto>>> ComplexBuffers;

public:
	Alignment(boost::shared_ptr<clArgStore> clArgs, boost::shared_ptr<boost::mutex> dlg_mtx) : clArguments(clArgs), _mtx(dlg_mtx) {};
	
	void StartAlign();
	void Process();
	void update_clArgs(boost::shared_ptr<clArgStore> clArgs) { ComplexBuffers.clear(); clArguments = clArgs; }
};