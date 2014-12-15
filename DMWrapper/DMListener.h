#include "DMListenable.h"
#include "DMROI.h"
#include "DMout.h"
#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>

#include "boost/thread/mutex.hpp"
#include "boost/thread.hpp"

// Templated so this can be called from any dialog box. Type still needs to have a Process method.

class DMListener : public boost::enable_shared_from_this<DMListener>
{
public:
	typedef DMListener this_type;
private:
	bool Active;
	std::string fPrefix;
	std::vector<DMListenable*> Listeners;
	std::vector<DMROI*> ROIs;

	boost::shared_ptr<DMListener> get_this_shared() { return shared_from_this(); }

	boost::mutex _mtx;
public:
	// so event handler can be closed
	ulong_ptr_t EventHandler;

	DMListener() : Active(false) {}
	DMListener(const char_range &prefix);

	~DMListener();

	// This is what will happen when the ROI changes.
	void OnROIChanged( uint32 img_disp_event_flags, const DM::ImageDisplay &img_disp, uint32 roi_change_flag, uint32 roi_disp_change_flags, const DM::ROI &roi );

	// This is what will happen when the data changes.
	void OnDataChanged( uint32 img_event_flags, const DM::Image &img);

	// Add a listener to a given image display registered to the OnROIChanged function.
	static boost::shared_ptr<DMListener> AddImageDisplayListener_Ptr(const DM::ImageDisplay &img_disp, const std::string &prefix);

	// Add a listener to a given image registered to the OnDataChanged function.
	static boost::shared_ptr<DMListener> AddImageListener_Ptr(const DM::Image &img, const std::string &prefix);

	void addListenable(DMListenable* listenable);

	void addROI(DMROI* roi);

	void AddDisplayEventListener(DM::ImageDisplay &disp);
	void AddImageEventListener(DM::Image &img);

	void RemoveDisplayEventListener(DM::ImageDisplay &disp);
	void RemoveImageEventListener(DM::Image &image);

	bool operator ==(const this_type &) const { return false; }

	bool IsActive() { return Active; }
	void SetActive(bool state) { Active = state; }

	void ClearListeners() { Listeners.clear(); }
	void ClearROIs() { ROIs.clear(); }
};