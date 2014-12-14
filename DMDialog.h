#pragma once

#include "stdafx.h"
#include "WorkerClass.h"
#include "clWrapper/clWrapper.h"
#include "clArgStore.h"
#include "Alignment.h"

#include <boost/shared_ptr.hpp>
#include "boost/thread/mutex.hpp"
#include "boost/thread.hpp"

class CDMDialog : public CDialog
{
	DECLARE_DYNCREATE(CDMDialog)

public:
	CDMDialog(CWnd* pParent = NULL);  // standard constructor
	virtual ~CDMDialog();

	enum { IDD = IDD_DMDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);  // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL Create(UINT templateID, CWnd* pParentWnd = NULL);
	
	virtual BOOL OnInitDialog();

	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnCbnSelchangeCmbDevices();
	afx_msg void OnBnClickedBtnGpa();
private:
	CComboBox combo_CLdev;
	boost::shared_ptr<boost::mutex> _mtx;
	std::list<clDevice> Devices;

	boost::shared_ptr<Alignment> Align;

	boost::shared_ptr<clArgStore> clArguments;
};
