#include "stdafx.h"

#include "Alignment.h"
#include<time.h>

#include"LinearAlgebra/Matrix.h"
#include"LinearAlgebra/Decomposition.h"


void Alignment::Process()
{
	boost::lock_guard<boost::mutex> lock(*_mtx);
	if (clArguments->haveDevice)
		Start();
	else
		DMresult << "Must select OpenCL device to do live GPA" << DMendl;
}


void Alignment::StartAlign()
{
	// Get the front image
	try
		{ Image.fromFront(); }
	catch (const std::invalid_argument& e)
		{ DMresult << "ERROR: " << e.what() << DMendl; return; }

	//Later need to add ability to select ROI (needs to be power of 2 etc)

	width = Image.getWidth();
	height = Image.getHeight();

	// Create OpenCL stuff
	clArguments->FFT.reset(new clFourier(*(clArguments->Context), width, height));
	if (!clArguments->CheckStatus("Setting up FFT")) return;
	for (int i = 0; i < 4; i++)
	{
		ComplexBuffers.push_back((*(clArguments->Context)).CreateBuffer<std::complex<float>, Auto>(width * height));
		if (!clArguments->CheckStatus("Creating buffer")) return;
	}

	Process();
}

void Alignment::DoWork()
{
	boost::lock_guard<boost::mutex> lock(WorkLock);

	width = Image.getWidth();
	height = Image.getHeight();
	depth = Image.getDepth();

	if (depth < 2)
	{
		DMresult << "ERROR: Image must be a stack." << DMendl;
		return;
	}

	OverDeterminedAlign();
}

void Alignment::OverDeterminedAlign()
{
	std::vector<std::complex<float>> data1(width*height);
	std::vector<std::complex<float>> data2(width*height);

	int n = depth;
	int m = (n - 1) * n / 2; // check this hasn't got a stupid rounding thing going on
	std::vector<float> b_x(m);
	std::vector<float> b_y(m);
	std::vector<float> r_x(m);
	std::vector<float> r_y(m);
	Matrix<float> A(m, n-1);
	int ind = 0;

	clock_t startTime = clock();

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			// do first cross-correlation and get pixel shift
			Image.GetData(data1, 0, 0, height, width, i, i + 1);
			Image.GetData(data2, 0, 0, height, width, j, j + 1);
			std::vector<std::complex<float>> XCF = CrossCorrelation(data1, data2);
			coord<int> pix = FindMaxima(XCF);

			// Refine the position
			std::vector<float> area(9);
			int max_ind = pix.y*width + pix.x;
			int ind2 = 0;
			for (int k = -1; k <= 1; k++)
				for (int m = -1; m <= 1; m++)
				{
					area[ind2] = XCF[(max_ind + k*width) + m].real();
					ind2++;
				}
			coord<float> sub_pix = FindVertexParabola(area);


			b_x[ind] = (static_cast<float>(pix.x) + sub_pix.x) - static_cast<float>(width) / 2;
			b_y[ind] = static_cast<float>(height) / 2 - (static_cast<float>(pix.y) + sub_pix.y);

			if (j == i + 1)
			{
				r_x[i] = b_x[ind];
				r_y[i] = b_y[ind];
			}

			// Create the coefficient matric here
			for (int k = i; k < j; k++)
				A(ind, k) = 1.0;


			ind++;
		}
	}

	std::vector<float> s_x = lsSolver(A, b_x);
	std::vector<float> s_y = lsSolver(A, b_y);

	for (int i = 0; i < s_x.size(); i++)
		DMresult << "(" << r_x[i] << ", " << r_y[i] << ") -> (" << s_x[i] << ", " << s_y[i] << ")" << DMendl;

	AlignImage(s_x, s_y);
}

std::vector<std::complex<float>> Alignment::CrossCorrelation(std::vector<std::complex<float>> image1, std::vector<std::complex<float>>image2)
{
	clWorkGroup GlobalWork(width, height, 1);

	ComplexBuffers[0]->Write(image1);
	ComplexBuffers[1]->Write(image2);

	(*(clArguments->FFT))(ComplexBuffers[0], ComplexBuffers[0], Direction::Forwards);
	(*(clArguments->FFT))(ComplexBuffers[1], ComplexBuffers[1], Direction::Forwards);

	clArguments->kMultiCorrelation->SetArg(0, ComplexBuffers[0], ArgumentType::Input);
	clArguments->kMultiCorrelation->SetArg(1, ComplexBuffers[1], ArgumentType::Input);
	clArguments->kMultiCorrelation->SetArg(2, ComplexBuffers[2], ArgumentType::Output);
	clArguments->kMultiCorrelation->SetArg(3, ComplexBuffers[3], ArgumentType::Output);
	clArguments->kMultiCorrelation->SetArg(4, width);
	clArguments->kMultiCorrelation->SetArg(5, height);

	(*clArguments->kMultiCorrelation)(GlobalWork);

	(*(clArguments->FFT))(ComplexBuffers[2], ComplexBuffers[2], Direction::Inverse);
	(*(clArguments->FFT))(ComplexBuffers[3], ComplexBuffers[3], Direction::Inverse);

	clArguments->kFFTShift->SetArg(0, ComplexBuffers[3], ArgumentType::Input);
	clArguments->kFFTShift->SetArg(1, ComplexBuffers[1], ArgumentType::Output);
	clArguments->kFFTShift->SetArg(2, width);
	clArguments->kFFTShift->SetArg(3, height);

	(*clArguments->kFFTShift)(GlobalWork);

	clArguments->Context->WaitForQueueFinish();

	//std::vector<std::complex<float>> PCF = ComplexBuffers[0]->GetLocal();
	return ComplexBuffers[1]->GetLocal();
}

coord<int> Alignment::FindMaxima(std::vector<std::complex<float>> data)
{
	float maxheight = 0.0f;
	int maxPosition = 0;

	for (int j = 0; j< width*height; j++)
	{
		// Use real part of ifft to get peak heights, not absolute....
		float val = data[j].real();

		if (val > maxheight)
		{
			maxheight = val;
			maxPosition = j;
		}
	}

	int y = maxPosition / width; // should be rounded down
	int x = maxPosition % (y * width);

	//return coord<int>((height / 2) - y, x - (width / 2));
	return coord<int>(y, x);
}

coord<float> Alignment::FindVertexParabola(std::vector<float> data)
{
	assert(data.size() == 9);

	float x1 = -1;
	float x2 = 0;
	float x3 = 1;

	float y1 = data[3];
	float y2 = data[4];
	float y3 = data[5];


	float denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	float A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	float B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	float C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	float xoffset = std::max(std::min(-B / (2 * A), 1.0f), -1.0f); // Bad fits cause subshifts that are way bigger than +-1.0
	//float peak1 = C - B*B / (4 * A);

	// centre/x values are the same
	y1 = data[1];
	y3 = data[7];

	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	float yoffset = std::max(std::min(-B / (2 * A), 1.0f), -1.0f);
	//float peak2 = C - B*B / (4 * A);

	//float peakheight = (peak1 + peak2) / (2.0f);

	return coord<float>(xoffset, yoffset);
}

void Alignment::AlignImage(std::vector<float> shiftx, std::vector<float> shifty)
{
	std::vector<float> cumulative_x(shiftx.size());
	std::vector<float> cumulative_y(shifty.size());

	cumulative_x[0] = shiftx[0];
	cumulative_y[0] = shifty[0];

	for (int i = 1; i < shiftx.size(); i++)
	{
		cumulative_x[i] = cumulative_x[i - 1] + shiftx[i];
		cumulative_y[i] = cumulative_y[i - 1] + shifty[i];
	}

	clWorkGroup GlobalWork(width, height, 1);
	std::vector<std::complex<float>> data(width*height);

	std::vector<float> summed(width*height);
	Image.GetData(summed, 0, 0, height, width, 0, 1);

	std::vector<DMImage> test_images;

	for (int i = 1; i < depth; i++)
	{
		Image.GetData(data, 0, 0, height, width, i, i + 1);
		ComplexBuffers[0]->Write(data);

		/*(*(clArguments->FFT))(ComplexBuffers[0], ComplexBuffers[0], Direction::Forwards);

		clArguments->kFFTShift->SetArg(0, ComplexBuffers[0], ArgumentType::Input);
		clArguments->kFFTShift->SetArg(1, ComplexBuffers[1], ArgumentType::Output);
		clArguments->kFFTShift->SetArg(2, width);
		clArguments->kFFTShift->SetArg(3, height);

		(*clArguments->kFFTShift)(GlobalWork);*/

		//clArguments->kShiftImage->SetArg(0, ComplexBuffers[1], ArgumentType::InputOutput);
		//clArguments->kShiftImage->SetArg(1, cumulative_x[i - 1]);
		//clArguments->kShiftImage->SetArg(2, cumulative_y[i - 1]);
		//clArguments->kShiftImage->SetArg(3, width);
		//clArguments->kShiftImage->SetArg(4, height);

		//(*clArguments->kShiftImage)(GlobalWork);

		clArguments->kBilinearInterpolate->SetArg(0, ComplexBuffers[0], ArgumentType::Input);
		clArguments->kBilinearInterpolate->SetArg(1, ComplexBuffers[1], ArgumentType::Output);
		clArguments->kBilinearInterpolate->SetArg(2, width);
		clArguments->kBilinearInterpolate->SetArg(3, height);
		clArguments->kBilinearInterpolate->SetArg(4, 0);
		clArguments->kBilinearInterpolate->SetArg(5, 0);
		clArguments->kBilinearInterpolate->SetArg(6, 0);
		clArguments->kBilinearInterpolate->SetArg(7, 0);
		clArguments->kBilinearInterpolate->SetArg(8, -cumulative_x[i - 1]);
		clArguments->kBilinearInterpolate->SetArg(9, -cumulative_y[i - 1]);
		clArguments->kBilinearInterpolate->SetArg(10, width);
		clArguments->kBilinearInterpolate->SetArg(11, height);
		clArguments->kBilinearInterpolate->SetArg(12, 0);
		clArguments->kBilinearInterpolate->SetArg(13, 0);

		(*clArguments->kBilinearInterpolate)(GlobalWork);


		/*clArguments->kFFTShift->SetArg(0, ComplexBuffers[1], ArgumentType::Input);
		clArguments->kFFTShift->SetArg(1, ComplexBuffers[0], ArgumentType::Output);

		(*clArguments->kFFTShift)(GlobalWork);

		(*(clArguments->FFT))(ComplexBuffers[0], ComplexBuffers[0], Direction::Inverse);*/

		std::vector<std::complex<float>> temp = ComplexBuffers[1]->GetLocal();

		//std::stringstream ss;
		//ss << "aligned " << i;
		//
		//test_images.push_back( DMImage(temp, ss.str().c_str(), 4, width, height) );

		for (int j = 0; j < temp.size(); j++)
			summed[j] += temp[j].real();
	}

	DMImage summ_image(summed, "summed image", 4, width, height);
}