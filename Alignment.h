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

class CDMDialog;

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

	// Reference to parent dialog. Used to update progress bar
	CDMDialog* parent;

	// Image used as input
	DMImage Image;

	// Image used to store input with blank frames removed
	DMImage BlankCorrected;

	// Dimensions of the dataset
	int width, height, depth;

	// Factor to use when applying exponential low pass filter
	float Bfactor;

	// Threshold for pixel movement when testing for good shifts
	float threshold;

	// Class storing all the needed OpenCL arguments (context, kernels etc.)
	boost::shared_ptr<clArgStore> clArguments;

	// Vector of complex float buffers used for OpenCL. Same size as input image.
	std::vector<boost::shared_ptr<clMemory<std::complex<float>, Auto>>> ComplexBuffers;

	// Complex flaot buffer used to get cropped output image. Different dimensions to image.
	boost::shared_ptr<clMemory<std::complex<float>, Auto>> OutputBuffer;

public:
	Alignment(){}
	Alignment(CDMDialog* prnt, boost::shared_ptr<clArgStore> clArgs, boost::shared_ptr<boost::mutex> dlg_mtx, float Bf, float thresh) :
		parent(prnt), clArguments(clArgs), _mtx(dlg_mtx), Bfactor(Bf), threshold(thresh)
	{
		DMresult << "ima doing da constructin" << DMendl;
	}

	//~Alignment()
	//{
	//	DMresult << "Align destroy 1" << DMendl;
	//	parent = 0;
	//	DMresult << "Align destroy 2" << DMendl;
	//	ComplexBuffers.clear();
	//	DMresult << "Align destroy 3" << DMendl;
	//	clArguments.reset();
	//	DMresult << "Align destroy 4" << DMendl;
	//	OutputBuffer.reset();
	//	DMresult << "Align destroy 5" << DMendl;
	//	_mtx.reset();
	//	DMresult << "Align destroy 6" << DMendl;
	//	Image.~DMImage();
	//	DMresult << "Align destroy 7" << DMendl;
	//	BlankCorrected.~DMImage();
	//	DMresult << "Align destroy 8" << DMendl;
	//}

	void StartAlign();
	void Process();
};