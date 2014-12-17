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

#include"LinearAlgebra/Matrix.h"

class Alignment : public WorkerClass, public DMListenable
{
private:
	// Main body of the alignment code. Runs in new thread to avoid locking anything else.
	void DoWork();

	//
	void RemoveBlankFrames();

	// Performs an alignment by solving an overdetermined set of equations relating all shifts
	void OverDeterminedAlign();

	// Performs Cross correlation and phase correlation using OpenCL
	std::vector<std::complex<float>> CrossCorrelation(std::vector<std::complex<float>> image1, std::vector<std::complex<float>>image2);

	// Finds the maxima of a vector and returns the x,y, position
	coord<int> FindMaxima(std::vector<std::complex<float>> &data);

	// Finds the minimum of a vector and returns the index
	int FindMinimaIndex(std::vector<float> &data);

	// Fit parabola to find refined maxima
	coord<float>FindVertexParabola(std::vector<float> data);

	// Finds values with large error and removes them from calculation
	std::vector<int> OverDeterminedThreshold(Matrix<std::complex<float>> &A, std::vector<std::complex<float>> &b, std::vector<float> &error, float thresh);

	// Method that performs the shifts and summation of the image
	void AlignImage(std::vector<float> shiftx, std::vector<float> shifty);

	// Mutex shared with the dialog class
	boost::shared_ptr<boost::mutex> _mtx;
	// Mutex purely for the Alignment class
	boost::mutex WorkLock;

	// Image used as input
	DMImage Image;

	// Image used to store input with blank frames removed
	DMImage BlankCorrected;

	// Dimensions of the dataset
	int width, height, depth;

	// Class storing all the needed OpenCL arguments (context, kernels etc.)
	boost::shared_ptr<clArgStore> clArguments;

	// Vector of complex float buffers used for OpenCL. Same size as input image.
	std::vector<boost::shared_ptr<clMemory<std::complex<float>, Auto>>> ComplexBuffers;

	// Complex flaot buffer used to get cropped output image. Different dimensions to image.
	boost::shared_ptr<clMemory<std::complex<float>, Auto>> OutputBuffer;

public:
	Alignment(){}
	Alignment(boost::shared_ptr<clArgStore> clArgs, boost::shared_ptr<boost::mutex> dlg_mtx) : clArguments(clArgs), _mtx(dlg_mtx) {};
	
	void StartAlign();
	void Process();
};