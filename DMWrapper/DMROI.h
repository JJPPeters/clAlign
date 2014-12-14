#pragma once
#include "../WorkerClass.h"
#include "DMListenable.h"

#include "boost/thread/mutex.hpp"
#include "boost/thread.hpp"

enum SCALE_INT{
	SCALE_ALL,
	SCALE_CENTRE,
	SCALE_TC,
	SCALE_TL,
	SCALE_LC,
	SCALE_BL,
	SCALE_BC,
	SCALE_BR,
	SCALE_RC,
	SCALE_TR
};

template <class T> class coord
{
public:
	T x, y;
	coord() : x(0), y(0) {}
	coord(T xx, T yy) : x(xx), y(yy) {}
};

template <class T> class rect
{
public:
	boost::shared_ptr<boost::mutex> _mtx;

	T t, l, b, r;
	rect() : t(0), l(0), b(0), r(0), _mtx(new boost::mutex) {}
	rect(T tt, T ll, T bb, T rr) : t(tt), l(ll), b(bb), r(rr), _mtx(new boost::mutex) {}
	rect(const rect<T> &copy) : t(copy.t), l(copy.l), b(copy.b), r(copy.r), _mtx(new boost::mutex) {}

	template <typename U> coord<U> centre() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>((l + r) / 2, (t + b) / 2); }
	template <typename U> coord<U> tc() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>((l + r) / 2, t); }
	template <typename U> coord<U> tl() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(l, t); }
	template <typename U> coord<U> lc() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(l, (t + b) / 2); }
	template <typename U> coord<U> bl() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(l, b); }
	template <typename U> coord<U> bc() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>((l + r) / 2, t); }
	template <typename U> coord<U> br() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(r, b); }
	template <typename U> coord<U> rc() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(r, (t + b) / 2); }
	template <typename U> coord<U> tr() { boost::lock_guard<boost::mutex> lock(*_mtx); return coord<U>(r, t); }

	T sizeX() { boost::lock_guard<boost::mutex> lock(*_mtx); return r - l; }
	T sizeY() { boost::lock_guard<boost::mutex> lock(*_mtx); return b - t; }
	T Area() { return sizeX() * sizeY(); }

	template <class U> bool compareSize(rect<U>& arg) { return sizeX() == arg.sizeX() && sizeY() == arg.sizeY(); }
	template <class U> bool operator==(rect<U>& RHS) { boost::lock_guard<boost::mutex> lock(*_mtx); return (t == RHS.t && b == RHS.b && l == RHS.l && r == RHS.r); }
	template <class U> void operator=(rect<U>& RHS) { boost::lock_guard<boost::mutex> lock(*_mtx); t = RHS.t; b = RHS.b; l = RHS.l; r = RHS.r;; }
};

struct DMROI: public WorkerClass, public DMListenable, public rect<int>
{
	boost::mutex process_mtx;

	//Holds the ROI
	DM::ROI roi;

	bool isPow2();

	void setRect();
	void setRect(int top, int left, int right, int bottom);
	void setRect(float top, float left, float right, float bottom);

	void scaleRect(float c);
	void scaleRect(float c, SCALE_INT choice);

	void updateRect();

	void addROI(DM::ImageDisplay &ImDisp);

	void fromImage(DM::ImageDisplay &ImDisp);
	void fromNew(DM::ImageDisplay &ImDisp, int top, int left, int bottom, int right, bool deletable, bool resizable);

	DMROI() {_mtx = boost::shared_ptr<boost::mutex>(new boost::mutex);}

	DMROI(DM::ImageDisplay &ImDisp);

	DMROI(DM::ImageDisplay &ImDisp, int top, int left, int bottom, int right, bool deletable, bool resizable);
	
	DMROI(const DMROI &copy) : roi(copy.roi) { _mtx = boost::shared_ptr<boost::mutex>(new boost::mutex); t = copy.t; l = copy.l; b = copy.b; r = copy.r; }

	~DMROI(){}

	void Process();

private:
	void DoWork();
};