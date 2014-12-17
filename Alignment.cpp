#include "stdafx.h"

#include "Alignment.h"
#include<time.h>
#include <numeric> // for accumulate function
#include <complex>

#include"LinearAlgebra/Decomposition.h"
#include"LinearAlgebra/SVD.h"


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

	// Should probably find a way of skipping this step if not needed?
	RemoveBlankFrames();

	OverDeterminedAlign();
}

void Alignment::RemoveBlankFrames()
{
	std::vector<float> data(width*height);
	std::vector<int> validZ(depth);
	int newDepth = 0;

	for (int i = 0; i < depth; i++)
	{
		Image.GetData(data, 0, 0, height, width, i, i + 1);

		if (data[0] != 0 || data[width * height - 1] != 0) // Quick test that should cover most non-zero cases
		{
			validZ[i] = 1;
			newDepth++;
		}
		else
		{
			// check rest of elements are also 0
			// might be somewhat intensive work?
			float sum_of_elems = std::accumulate(data.begin(), data.end(), 0);

			if (sum_of_elems == 0)
				validZ[i] = 0;
			else
			{
				validZ[i] = 1;
				newDepth++;
			}
		}
	}

	// Create the new image now.
	BlankCorrected = DMImage("Blank frame corrected", 4, width, height, newDepth);

	int ind = 0;

	for (int i = 0; i < depth; i++)
	{
		if (validZ[i] == 1)
		{
			// get data then set it.
			try
				{ Image.GetData(data, 0, 0, height, width, i, i + 1); }
			catch (const std::invalid_argument& e)
				{ DMresult << "ERROR: " << e.what() << DMendl; break; }
			try
				{ BlankCorrected.SetRealData(data, width*height*ind); }
			catch (const std::invalid_argument& e)
				{ DMresult << "ERROR: " << e.what() << DMendl; break; }
			ind++;
		}
	}

	// update depth
	depth = newDepth;
}

void Alignment::OverDeterminedAlign()
{
	std::vector<std::complex<float>> data1(width*height);
	std::vector<std::complex<float>> data2(width*height);

	int n = depth;
	int m = (n - 1) * n / 2; // check this hasn't got a stupid rounding thing going on
	std::vector<std::complex<float>> b(m);
	std::vector<std::complex<float>> r(n-1);
	Matrix<std::complex<float>> A(m, n - 1);
	int ind = 0;

	clock_t startTime = clock();

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			// do first cross-correlation and get pixel shift
			BlankCorrected.GetData(data1, 0, 0, height, width, i, i + 1);
			BlankCorrected.GetData(data2, 0, 0, height, width, j, j + 1);

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

			// reset origin and apply refinement
			float t1 = (static_cast<float>(pix.x) + sub_pix.x) - static_cast<float>(width) / 2;
			float t2 = static_cast<float>(height) / 2 - (static_cast<float>(pix.y) + sub_pix.y);

			b[ind] = std::complex<float>(t1, t2);

			if (j == i + 1)
				r[i] = b[ind];

			// Create the coefficient matrix here
			for (int k = i; k < j; k++)
				A(ind, k) = 1.0;

			ind++;
		}
	}

	std::vector<std::complex<float>> s = lsSolver(A, b);

	std::vector<float> s_x(m);
	std::vector<float> s_y(m);

	for (int i = 0; i < s.size(); i++)
	{
		s_x[i] = s[i].real();
		s_y[i] = s[i].imag();
	}

	std::vector<std::complex<float>> temp = A*s;

	std::vector<float> error(temp.size());

	for (int i = 0; i < temp.size(); ++i)
		error[i] = std::abs(temp[i] - b[i]);

	std::vector<int> goodlist = OverDeterminedThreshold(A, b, error, 5);

	std::vector<float> s2_x(m);
	std::vector<float> s2_y(m);

	s = lsSolver(A, b);

	for (int i = 0; i < s.size(); i++)
	{
		DMresult << "(" << r[i].real() << ", " << r[i].imag() << ") -> (" << s_x[i] << ", " << s_y[i] << ") -> (" << s[i].real() << ", " << s[i].imag() << ")" << DMendl;
		s2_x[i] = s[i].real();
		s2_y[i] = s[i].imag();
	}

	AlignImage(s2_x, s2_y);
}

std::vector<std::complex<float>> Alignment::CrossCorrelation(std::vector<std::complex<float>> image1, std::vector<std::complex<float>>image2)
{
	clWorkGroup GlobalWork(width, height, 1);

	ComplexBuffers[0]->Write(image1);
	ComplexBuffers[1]->Write(image2);

	(*(clArguments->FFT))(ComplexBuffers[0], ComplexBuffers[0], Direction::Forwards);
	(*(clArguments->FFT))(ComplexBuffers[1], ComplexBuffers[1], Direction::Forwards);

	// Needs to be user set at some point
	float B = 100.0;

	clArguments->kExponentialPass->SetArg(0, ComplexBuffers[0], ArgumentType::InputOutput);
	clArguments->kExponentialPass->SetArg(1, B);
	clArguments->kExponentialPass->SetArg(2, width);
	clArguments->kExponentialPass->SetArg(3, height);

	(*clArguments->kExponentialPass)(GlobalWork);

	clArguments->kExponentialPass->SetArg(0, ComplexBuffers[1], ArgumentType::InputOutput);

	(*clArguments->kExponentialPass)(GlobalWork);

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

	return ComplexBuffers[1]->GetLocal();
}

coord<int> Alignment::FindMaxima(std::vector<std::complex<float>> &data)
{
	assert(data.size() = width * height);

	float maxValue = 0.0f;
	int maxPosition = 0;

	for (int j = 0; j< width*height; j++)
	{
		// Use real part of ifft to get peak heights, not absolute....
		float val = data[j].real();

		if (val > maxValue)
		{
			maxValue = val;
			maxPosition = j;
		}
	}

	int y = maxPosition / width; // should be rounded down
	int x = maxPosition % (y * width);

	return coord<int>(x, y);
}

int Alignment::FindMinimaIndex(std::vector<float> &data)
{
	int minPosition = 0;
	float minValue = data[0];
	for (int i = 1; i<data.size(); i++)
	{
		if (data[i]<minValue)
		{
			minValue = data[i];
			minPosition = i;
		}
	}

	return minPosition;
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

	if (denom == 0)
		return coord<float>(0.0, 0.0);

	float A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	float B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	float C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	if (A == 0)
		return coord<float>(0.0, 0.0);

	float xoffset = std::max(std::min(-B / (2 * A), 1.0f), -1.0f); // Bad fits cause subshifts that are way bigger than +-1.0

	y1 = data[1];
	y3 = data[7];

	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

	if (denom == 0)
		return coord<float>(0.0, 0.0);

	A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	if (A == 0)
		return coord<float>(0.0, 0.0);

	float yoffset = std::max(std::min(-B / (2 * A), 1.0f), -1.0f);

	return coord<float>(xoffset, yoffset);
}

std::vector<int> Alignment::OverDeterminedThreshold(Matrix<std::complex<float>> &A, std::vector<std::complex<float>> &b, std::vector<float> &error, float thresh)
{
	std::vector<int> goodlist, bk;
	Matrix<std::complex<float>> tA = A;
	std::vector<std::complex<float>> tb = b;
	std::vector<float> terror = error;
	SVD<std::complex<float>> svd;

	Matrix<float> temp(10,10);

	int rank = 0;
	int i, pos, ii;
	for (i = 0; i<terror.size(); i++)
	{
		pos = FindMinimaIndex(terror);
		goodlist.push_back(pos);
		terror[pos] = 1e20;
	}
	bk = goodlist;

	for (i = goodlist.size() - 1; i >= A.cols(); i--)
	{
		if (error[goodlist[i]]<thresh)
			break;
	}
	ii = i;

	do
	{
		goodlist = bk;
		goodlist.erase(goodlist.begin() + ii + 1, goodlist.end());
		sort(goodlist.begin(), goodlist.end());

		A = Matrix<std::complex<float> >(goodlist.size(), A.cols(), 0.0);
		b = std::vector<std::complex<float> >(goodlist.size(), 0.0);
		for (i = 0; i<goodlist.size(); i++)
		{
			A.setRow(tA.getRow(goodlist[i]), i);
			b[i] = tb[goodlist[i]];
		}

		svd.decompose(A);
		rank = svd.rank();

		ii++;
	} while (rank<A.cols() && ii<bk.size());

	return goodlist;
}

void Alignment::AlignImage(std::vector<float> shiftx, std::vector<float> shifty)
{
	// First calculate the cumulative shifts

	// Vectors to contain cumulative shifts
	std::vector<float> cumulative_x(shiftx.size()+1);
	std::vector<float> cumulative_y(shifty.size()+1);

	// Stores maximum shifts so we can crop the images appropriately
	float max_x = 0;
	float max_y = 0;
	float min_x = 0;
	float min_y = 0;

	// First shift can be set manually
	cumulative_x[0] = 0;
	cumulative_y[0] = 0;

	// Loop through and add shifts
	for (int i = 1; i < shiftx.size()+1; i++)
	{
		cumulative_x[i] = cumulative_x[i-1] + shiftx[i-1];
		cumulative_y[i] = cumulative_y[i-1] + shifty[i-1];

		// Test for new min/max
		if (cumulative_x[i] > max_x)
			max_x = cumulative_x[i];
		if (cumulative_y[i] > max_y)
			max_y = cumulative_y[i];
		if (cumulative_x[i] < min_x)
			min_x = cumulative_x[i];
		if (cumulative_y[i] < min_y)
			min_y = cumulative_y[i];
	}

	// Convert max/min shifts to integers
	int iMax_x = static_cast<int>(floor(max_x));
	int iMax_y = static_cast<int>(floor(max_y));
	int iMin_x = static_cast<int>(floor(min_x));
	int iMin_y = static_cast<int>(floor(min_y));

	// Calculate new image sizes from drifts
	int newWidth = width - (iMax_x - iMin_x);
	int newHeight = height - (iMax_y - iMin_y);

	OutputBuffer = (*(clArguments->Context)).CreateBuffer<std::complex<float>, Auto>(newWidth * newHeight);

	clWorkGroup GlobalWork(width, height, 1);
	std::vector<std::complex<float>> data(width*height);

	std::vector<float> summed(newWidth*newHeight);

	for (int i = 0; i < newWidth*newHeight; i++)
		summed[i] = 0.0;

	DMImage aligned_image("Aligned Image", 4, newWidth, newHeight, depth);

	for (int i = 0; i < depth; i++)
	{
		BlankCorrected.GetData(data, 0, 0, height, width, i, i + 1);
		ComplexBuffers[0]->Write(data);

		// Shift image with no padding
		clArguments->kBilinearInterpolate->SetArg(0, ComplexBuffers[0], ArgumentType::Input);
		clArguments->kBilinearInterpolate->SetArg(1, OutputBuffer, ArgumentType::Output);
		clArguments->kBilinearInterpolate->SetArg(2, width);
		clArguments->kBilinearInterpolate->SetArg(3, height);
		clArguments->kBilinearInterpolate->SetArg(4, 0);
		clArguments->kBilinearInterpolate->SetArg(5, 0);
		clArguments->kBilinearInterpolate->SetArg(6, 0);
		clArguments->kBilinearInterpolate->SetArg(7, 0);
		clArguments->kBilinearInterpolate->SetArg(8, cumulative_x[i]);
		clArguments->kBilinearInterpolate->SetArg(9, -cumulative_y[i]);
		clArguments->kBilinearInterpolate->SetArg(10, newWidth);
		clArguments->kBilinearInterpolate->SetArg(11, newHeight);

		//int temp1 = iMax_y - static_cast<int>(floor(cumulative_y[i]));
		//int temp2 = iMax_x - static_cast<int>(floor(cumulative_x[i]));
		int temp1 = 0;
		int temp2 = 0;

		clArguments->kBilinearInterpolate->SetArg(12, temp1 );
		clArguments->kBilinearInterpolate->SetArg(13, temp2 );

		(*clArguments->kBilinearInterpolate)(GlobalWork);

		std::vector<std::complex<float>> temp = OutputBuffer->GetLocal();

		try
		{
			aligned_image.SetComplexData(temp, newWidth*newHeight*i, SHOW_REAL);
		}
		catch (const std::invalid_argument& e)
		{
			DMresult << "ERROR: " << e.what() << DMendl; break;
		}

		for (int j = 0; j < summed.size(); j++)
			summed[j] = summed[j] + temp[j].real();

	}

	DMImage summ_image(summed, "Summed image", 4, newWidth, newHeight);
}