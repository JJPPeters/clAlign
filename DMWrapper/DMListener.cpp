#include "stdafx.h"
#include "DMListener.h"

DMListener::DMListener(const char_range &prefix) : fPrefix(prefix.begin(), prefix.end()), Active(true) {}

DMListener::~DMListener(){}

void DMListener::OnROIChanged(uint32 img_disp_event_flags, const DM::ImageDisplay &img_disp, uint32 roi_change_flag, uint32 roi_disp_change_flags, const DM::ROI &roi)
{
	// BE CAREFUL!! The listener spams out 3 or 4 signals each time. This should be safe as the Process method checks to see if it is already working.
	for (int i = 0; i < ROIs.size(); i++)
	{
		if ((*(ROIs[i])).roi.GetID() == roi.GetID())
		{
			boost::lock_guard<boost::mutex> lock(_mtx);
			(ROIs[i])->Process();
		}
	}

	for (int i = 0; i < Listeners.size(); i++)
	{
		boost::lock_guard<boost::mutex> lock(_mtx);
		(Listeners[i])->Process();
	}
}

void DMListener::OnDataChanged(uint32 img_event_flags, const DM::Image &img)
{
	// BE CAREFUL!! The listener spams out 3 or 4 signals each time. This should be safe as the Process method checks to see if it is already working.
	for (int i = 0; i < Listeners.size(); i++)
	{
		boost::lock_guard<boost::mutex> lock(_mtx);
		(Listeners[i])->Process();
	}
}

boost::shared_ptr<DMListener> DMListener::AddImageDisplayListener_Ptr(const DM::ImageDisplay &img_disp, const std::string &prefix)
{
	boost::shared_ptr<this_type> ptr(new this_type(prefix));

	ptr->EventHandler = DM::ImageDisplayAddEventListener
		(
		img_disp
		, DM::NewAnonymousCppPtrScriptObject(ptr
		, dmscript::function_list
		(
		dmscript::dms_bindfn_2<cpp_method_traits_5<this_type, void, uint32, const DM::ImageDisplay &, uint32, uint32, const DM::ROI &>
		, &this_type::OnROIChanged>("OnROIChanged", dmscript::dms_func_attr_threadsafe)
		)
		)
		, "roi_changed:OnROIChanged"
		, "liveGPA"
		);

	return ptr;
}

boost::shared_ptr<DMListener> DMListener::AddImageListener_Ptr(const DM::Image &img, const std::string &prefix)
{
	boost::shared_ptr<this_type> ptr(new this_type(prefix));

	ptr->EventHandler = DM::ImageAddEventListener
		(
		img
		, DM::NewAnonymousCppPtrScriptObject(ptr
		, dmscript::function_list
		(
		dmscript::dms_bindfn_2<cpp_method_traits_2<this_type, void, uint32, const DM::Image &>
		, &this_type::OnDataChanged>("OnDataChanged", dmscript::dms_func_attr_threadsafe)
		)
		)
		, "data_changed,slices_changed:OnDataChanged"
		, "liveGPA"
		);

	return ptr;
}

void DMListener::addListenable(DMListenable* listenable)
{
	Listeners.push_back(listenable);
}

void DMListener::addROI(DMROI* roi)
{
	ROIs.push_back(roi);
}

void DMListener::RemoveDisplayEventListener(DM::ImageDisplay &disp)
{
	if (Active)
	{
		DM::ImageDisplayRemoveEventListener(disp, EventHandler);
		Active = false;
	}
}

void DMListener::RemoveImageEventListener(DM::Image &image)
{
	if (Active)
	{
		DM::ImageRemoveEventListener(image, EventHandler);
		Active = false;
	}
}

void DMListener::AddDisplayEventListener(DM::ImageDisplay &disp)
{
	Active = false;

	boost::shared_ptr<this_type> ptr = get_this_shared();

	ptr->EventHandler = DM::ImageDisplayAddEventListener
		(
		disp
		, DM::NewAnonymousCppPtrScriptObject(ptr
		, dmscript::function_list
		(
		dmscript::dms_bindfn_2<cpp_method_traits_5<this_type, void, uint32, const DM::ImageDisplay &, uint32, uint32, const DM::ROI &>
		, &this_type::OnROIChanged>("OnROIChanged", dmscript::dms_func_attr_threadsafe)
		)
		)
		, "roi_changed:OnROIChanged"
		, "liveGPA"
		);

	Active = true;
}

void DMListener::AddImageEventListener(DM::Image &img)
{
	Active = false;

	boost::shared_ptr<this_type> ptr = get_this_shared();

	ptr->EventHandler = DM::ImageAddEventListener
		(
		img
		, DM::NewAnonymousCppPtrScriptObject(ptr
		, dmscript::function_list
		(
		dmscript::dms_bindfn_2<cpp_method_traits_2<this_type, void, uint32, const DM::Image &>
		, &this_type::OnDataChanged>("OnDataChanged", dmscript::dms_func_attr_threadsafe)
		)
		)
		, "data_changed,slices_changed:OnDataChanged"
		, "liveGPA"
		);

	Active = true;
}