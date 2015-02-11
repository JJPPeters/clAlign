#pragma once

#include "DMListener.h"
#include "DMROI.h"
#include "DMout.h"

#include <boost/shared_ptr.hpp>

enum ShowComplex{
	SHOW_REAL,
	SHOW_IMAG,
	SHOW_AMP,
	SHOW_PHASE,
	SHOW_ABS
};

class DMImage
{
private:
	unsigned long width;
	unsigned long height;
	unsigned long depth;
public:
	DM::Image Image;
	DM::ImageDisplay Display;

	bool isSquarePow2()
	{
		return ((width != 0) && !(width & (width - 1))) && ((height != 0) && !(height & (height - 1))) && (width == height);
	}

	unsigned long getWidth(){ return width; };
	unsigned long getHeight(){ return height; };
	unsigned long getDepth(){ return depth; };

	std::vector<boost::shared_ptr<DMROI>> ROIs;

	boost::shared_ptr<DMListener> ROIListener;
	boost::shared_ptr<DMListener> DataListener;

	DMImage() : ROIListener(new DMListener), DataListener(new DMListener) {}

	template <typename T>
	DMImage(const std::vector<std::complex<T>> &data, std::string title, int bytes, int x, int y, int z) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z);
		SetComplexData(data, SHOW_REAL);
	}

	template <typename T>
	DMImage(const std::vector<std::complex<T>> &data, std::string title, int bytes, int x, int y, int z, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z, t, l, b, r);
		SetComplexData(data, SHOW_REAL);
	}

	template <typename T>
	DMImage(const std::vector<T> &data, std::string title, int bytes, int x, int y, int z) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z);
		SetRealData(data);
	}

	template <typename T>
	DMImage(const std::vector<T> &data, std::string title, int bytes, int x, int y, int z, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z, t, l, b, r);
		SetRealData(data);
	}

	template <typename T>
	DMImage(const std::vector<std::complex<T>> &data, std::string title, int bytes, int x, int y) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y);
		SetComplexData(data, SHOW_REAL);
	}

	template <typename T>
	DMImage(const std::vector<std::complex<T>> &data, std::string title, int bytes, int x, int y, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, t, l, b, r);
		SetComplexData(data, SHOW_REAL);
	}

	template <typename T>
	DMImage(const std::vector<T> &data, std::string title, int bytes, int x, int y) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y);
		SetRealData(data);
	}

	template <typename T>
	DMImage(const std::vector<T> &data, std::string title, int bytes, int x, int y, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, t, l, b, r);
		SetRealData(data);
	}

	DMImage(std::string title, int bytes, int x, int y, int z) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z);
	}

	DMImage(std::string title, int bytes, int x, int y, int z, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, z, t, l, b, r);
	}

	DMImage(std::string title, int bytes, int x, int y) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y);
	}

	DMImage(std::string title, int bytes, int x, int y, long t, long l, long b, long r) : ROIListener(new DMListener), DataListener(new DMListener)
	{
		fromReal(title, bytes, x, y, t, l, b, r);
	}

	~DMImage() { RemoveDataListener(); RemoveROIListener(); }

	void fromFront();

	void fromReal(std::string title, int bytes, int x, int y, int z);
	void fromReal(std::string title, int bytes, int x, int y, int z, long t, long l, long b, long r);

	void fromReal(std::string title, int bytes, int x, int y);
	void fromReal(std::string title, int bytes, int x, int y, long t, long l, long b, long r);

	void Reshape(float Ratio, int x, int y);

	void DeleteImage() { DM::DeleteImage(Image); }

	void addROI(DMROI roi);

	void CreateROIListener();
	void AddROIListener();
	void RemoveROIListener();

	void CreateDataListener();
	void AddDataListener();
	void RemoveDataListener();

	void ClearListenables(){ (*ROIListener).ClearListeners(); (*ROIListener).ClearROIs(); }

	coord<long> GetWindowPosition() { long x, y; DM::GetWindowPosition(Image, &x, &y); return coord<long>(x, y); }
	coord<long> GetWindowSize() { long x, y; DM::GetWindowSize(Image, &x, &y); return coord<long>(x, y); }

	void SetWindowPosition(coord<long> pos) { DM::SetWindowPosition(Image, pos.x, pos.y); }
	void SetWindowSize(coord<long> size) { DM::SetWindowSize(Image, size.x, size.y); }

	template <typename T>
	void GetData(std::vector<T> &destination) { GetData(destination, 0, 0, height, width, 0, 1); }

	template <typename T>
	void GetData(std::vector<T> &destination, DMROI crop) { GetData(destination, crop.t, crop.l, crop.b, crop.r, 0, 1); }

	template <typename T>
	void GetData(std::vector<T> &destination, int t, int l, int b, int r, int front, int back);

	template <typename T>
	void SetRealData(const std::vector<T> &source, int offset);
	template <typename T>
	void SetRealData(const std::vector<T> &source);

	template <typename T>
	void SetComplexData(const std::vector<T> &source, int offset, ShowComplex doComplex);
	template <typename T>
	void SetComplexData(const std::vector<T> &source, ShowComplex doComplex);

private:

	template <typename U, typename T>
	void GetLockerData(std::vector<T> &destination, int t, int l, int b, int r, int front, int back);

	template <typename U, typename T>
	void SetRealLockerData(const std::vector<T> &source, int offset);

	template <typename U, typename T>
	void SetComplexLockerData(const std::vector<T> &source, int offeset, ShowComplex doComplex);

};

// This is where the template parts are implemented
#include "DMImage-ReadWrite.h"