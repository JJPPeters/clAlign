#include "stdafx.h"
#include "DialogPlugIn.h"
#include "DMPluginStubs.h"
#include "DMDialog.h"
#include "DMWrapper/DMout.h"

IMPLEMENT_DYNCREATE(CDMDialog, CDialog)

BEGIN_MESSAGE_MAP(CDMDialog, CDialog)
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_BN_CLICKED(IDC_BTN_GPA, &CDMDialog::OnBnClickedBtnGpa)
	ON_CBN_SELCHANGE(IDC_CMB_DEVICES, &CDMDialog::OnCbnSelchangeCmbDevices)
END_MESSAGE_MAP()

CDMDialog::CDMDialog(CWnd* pParent) : CDialog(CDMDialog::IDD, pParent), _mtx(new boost::mutex), clArguments(new clArgStore) {}

CDMDialog::~CDMDialog(){}

void CDMDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_CMB_DEVICES, combo_CLdev);
}

BOOL CDMDialog::Create(UINT templateID, CWnd* pParentWnd){return CDialog::Create(IDD, pParentWnd);}

void CDMDialog::OnPaint(){CDialog::OnPaint();}

BOOL CDMDialog::OnEraseBkgnd(CDC* pDC){return CDialog::OnEraseBkgnd(pDC);}

BOOL CDMDialog::OnInitDialog()
{
	CDialog::OnInitDialog();

	Devices = OpenCL::GetDeviceList();
	int numDevices = Devices.size();

	if (numDevices < 1)
		{ DMresult << "ERROR: " << "No OpenCL devices found" << DMendl; return false; }

	// Used to index position in list and is stored in combobox
	int ItemData = 0;

	for (std::list<clDevice>::iterator iterator = Devices.begin(), end = Devices.end(); iterator != end; ++iterator)
	{
		int i;
		clDevice dev = *iterator;
		i = combo_CLdev.AddString(dev.GetDeviceName().c_str());
		combo_CLdev.SetItemData(i, (DWORD)ItemData);
		ItemData++;
	}

	combo_CLdev.EnableWindow(true);

	return true;
}

void CDMDialog::OnCbnSelchangeCmbDevices()
{
	boost::lock_guard<boost::mutex> lock(*_mtx);

	DMresult << "Switching OpenCL device, this can take some time. Please wait..." << DMendl;

	int sel = combo_CLdev.GetCurSel();
	int item = combo_CLdev.GetItemData(sel);

	std::list<clDevice>::iterator iterator = Devices.begin();
	std::advance(iterator, item);
	clDevice dev = *iterator;

	DMresult << "Selected device: " << dev.GetDeviceName() << DMendl;

	clArguments.reset(new clArgStore(dev));

	DMresult << "Setting up complete, good to go!" << DMendl;
}

void CDMDialog::OnBnClickedBtnGpa()
{
	if (!(clArguments->haveDevice))
		{ DMresult << "ERROR: " << " No OpenCL device selected" << DMendl; return; }
	Align.reset(new Alignment(clArguments, _mtx));
	Align->StartAlign();
}
