#include "stdafx.h"

#include "DMImage.h"
#include "DMout.h"

void DMImage::fromFront()
{
	int success = DM::GetFrontImage(Image);

	if (!success)
		throw std::invalid_argument("Cannot get front image");
	// is there ever an index != 0?
	Display = Image.GetImageDisplay(0);
	Image.GetDimensionSizes(width, height, depth);
	Image.GetDimensionSizes();
	CreateDataListener();
}

void DMImage::fromReal(std::string title, int bytes, int x, int y, int z)
{
	Image = DM::RealImage(title, bytes, x, y, z);
	Image.GetOrCreateImageDocument().Show();
	Image.GetOrCreateImageDocument().Clean();
	Display = Image.GetImageDisplay(0);
	Image.GetDimensionSizes(width, height, depth);
	CreateDataListener();
}

void DMImage::fromReal(std::string title, int bytes, int x, int y, int z, long t, long l, long b, long r)
{
	Image = DM::RealImage(title, bytes, x, y, z);
	Image.GetOrCreateImageDocument().ShowAtRect(t, l, b, r);
	Image.GetOrCreateImageDocument().Clean();
	Display = Image.GetImageDisplay(0);
	Image.GetDimensionSizes(width, height, depth);
	CreateDataListener();
}

void DMImage::fromReal(std::string title, int bytes, int x, int y)
{
	fromReal(title, bytes, x, y, 1);
}

void DMImage::fromReal(std::string title, int bytes, int x, int y, long t, long l, long b, long r)
{
	fromReal(title, bytes, x, y, 1, t, l, b, r);
}

void DMImage::Reshape(float Ratio, int x, int y)
{
	if (depth > 1)
		DMresult << "WARNING: Reshaping not tested with stack images." << DMendl;

	coord<long> pos = GetWindowPosition();
	coord<long> size = GetWindowSize();
	
	RemoveDataListener(); // Causes crashes if listener hasn't been created
	RemoveROIListener();

	std::string Title = Image.GetName();
	long Bytes = Image.GetDataElementBitSize() / 8;

	DeleteImage();  // Doesnt delete ROIs but they need re-adding

	fromReal(Title, Bytes, x, y, pos.y, pos.x, pos.y + size.y, pos.x + size.x);

	AddROIListener();
	AddDataListener();

	for (int i = 0; i < ROIs.size(); i++)
	{
		ROIs[i]->scaleRect(Ratio);
		ROIs[i]->addROI(Display);
	}
}

void DMImage::addROI(DMROI roi)
{
	ROIs.push_back(boost::make_shared<DMROI>(roi));

	if (!(*ROIListener).IsActive())
		CreateROIListener();

	(*ROIListener).addROI(ROIs.back().get());
}

void DMImage::CreateROIListener()
{
	ROIListener = DMListener::AddImageDisplayListener_Ptr(Display, "ROIListener");
}

void DMImage::RemoveROIListener()
{
	ROIListener->RemoveDisplayEventListener(Display);
}

void DMImage::AddROIListener()
{
	ROIListener->AddDisplayEventListener(Display);
}

void DMImage::CreateDataListener()
{
	DataListener = DMListener::AddImageListener_Ptr(Image, "DataListener");
}

void DMImage::RemoveDataListener()
{
	DataListener->RemoveImageEventListener(Image);
}

void DMImage::AddDataListener()
{
	DataListener->AddImageEventListener(Image);
}