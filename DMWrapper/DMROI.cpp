#include "stdafx.h"
#include "DMROI.h"
#include "DMout.h"

void DMROI::Process(){ boost::lock_guard<boost::mutex> lock(process_mtx); Start(); }

void DMROI::DoWork(){ updateRect(); }

bool DMROI::isPow2()
{
	boost::lock_guard<boost::mutex> lock(*_mtx);
	return ((r-l != 0) && !(r-l & (r-l - 1))) && ((b-t != 0) && !(b-t & (b-t - 1))) && (r-l == b-t);
	//return ((sizeX() != 0) && !(sizeX() & (sizeX() - 1))) && ((sizeY() != 0) && !(sizeY() & (sizeY() - 1))) && (sizeX() == sizeY());
}

void DMROI::setRect()
{
	roi.SetRectangle(t, l, b, r);
}

void DMROI::setRect(int top, int left, int right, int bottom)
{
	t = top;
	l = left;
	r = right;
	b = bottom;

	setRect();
}

void DMROI::setRect(float top, float left, float right, float bottom)
{
	t = static_cast<int>(top);
	l = static_cast<int>(left);
	r = static_cast<int>(right);
	b = static_cast<int>(bottom);

	setRect();
}

void DMROI::scaleRect(float c)
{
	boost::lock_guard<boost::mutex> lock(*_mtx); 
	t = static_cast<int>(c*t);
	l = static_cast<int>(c*l);
	r = static_cast<int>(c*r);
	b = static_cast<int>(c*b);

	setRect();
}

void DMROI::scaleRect(float c, SCALE_INT choice)
{
	if (choice == SCALE_ALL){ scaleRect(c); return; }
	
	int newWidth = (int)(c * sizeX());
	int newHeight = (int)(c * sizeY());

	boost::lock_guard<boost::mutex> lock(*_mtx); 

	coord<int> origin;

	switch (choice)
	{
		case SCALE_CENTRE:
			origin = centre<int>(); break;
		case SCALE_TC:
			origin = tc<int>(); break;
		case SCALE_TL:
			origin = tl<int>(); break;
		case SCALE_LC:
			origin = lc<int>(); break;
		case SCALE_BL:
			origin = bl<int>(); break;
		case SCALE_BC:
			origin = bc<int>(); break;
		case SCALE_BR:
			origin = br<int>(); break;
		case SCALE_RC:
			origin = rc<int>(); break;
		case SCALE_TR:
			origin = tr<int>(); break;
	}

	t = c * (t - origin.y) + origin.y;
	b = c * (b - origin.y) + origin.y;

	l = c * (l - origin.x) + origin.x;
	r = c * (r - origin.x) + origin.x;

	setRect();
}

void DMROI::updateRect()
{
	boost::lock_guard<boost::mutex> lock(*_mtx);
	float top, left, right, bottom;
	roi.GetRectangle(&top, &left, &bottom, &right);
	t = static_cast<int>(top);
	l = static_cast<int>(left);
	b = static_cast<int>(bottom);
	r = static_cast<int>(right);
}

void DMROI::addROI(DM::ImageDisplay &ImDisp)
{
	ImDisp.AddROI(roi);
	roi.SetRectangle(t, l, b, r);
}

void DMROI::fromImage(DM::ImageDisplay &ImDisp)
{
	long num = ImDisp.CountROIs();
	if (num < 1)
		throw std::invalid_argument("Image contains no ROIs");

	int selectedID = -1;

	for (long i = 0; i < num; i++)
	{
		bool selected = ImDisp.IsROISelected(ImDisp.GetROI(i));
		if (selected)
			selectedID = i;
	}

	if (selectedID == -1)
		throw std::invalid_argument("No ROI in the image selected");

	roi = ImDisp.GetROI(selectedID);
	updateRect();
}

void DMROI::fromNew(DM::ImageDisplay &ImDisp, int top, int left, int bottom, int right, bool deletable, bool resizable)
{
	t = top;
	l = left;
	r = right;
	b = bottom;

	roi = DigitalMicrograph::NewROI();
	roi.SetRectangle(t, l, b, r);
	ImDisp.AddROI(roi);

	roi.SetDeletable(deletable);
	roi.SetResizable(resizable);
}

DMROI::DMROI(DM::ImageDisplay &ImDisp)
{
	 _mtx = boost::shared_ptr<boost::mutex>(new boost::mutex); 
	fromImage(ImDisp);
}

DMROI::DMROI(DM::ImageDisplay &ImDisp, int top, int left, int bottom, int right, bool deletable, bool resizable)
{
	 _mtx = boost::shared_ptr<boost::mutex>(new boost::mutex); 
	fromNew(ImDisp, top, left, bottom, right, deletable, resizable);
}