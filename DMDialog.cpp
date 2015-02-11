#include "stdafx.h"
#include "DialogPlugIn.h"
#include "DMPluginStubs.h"
#include "DMDialog.h"
#include "DMWrapper/DMout.h"

#include"boost/lexical_cast.hpp"

IMPLEMENT_DYNCREATE(CDMDialog, CDialog)

BEGIN_MESSAGE_MAP(CDMDialog, CDialog)
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_BN_CLICKED(IDC_BTN_ALIGN, &CDMDialog::OnBnClickedBtnGpa)
	//ON_BN_CLICKED(IDC_RDO_XCF, &CDMDialog::OnBnClickedBtnGpa)
	//ON_BN_CLICKED(IDC_RDO_PCF, &CDMDialog::OnBnClickedBtnGpa)
	ON_CBN_SELCHANGE(IDC_CMB_DEVICES, &CDMDialog::OnCbnSelchangeCmbDevices)
END_MESSAGE_MAP()

CDMDialog::CDMDialog(CWnd* pParent) : CDialog(CDMDialog::IDD, pParent), _mtx(new boost::mutex), s_bfactor("150"), s_thresh("5")
{}

CDMDialog::~CDMDialog(){}

void CDMDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_CMB_DEVICES, combo_CLdev);
	DDX_Control(pDX, IDC_CMB_METHOD, combo_method);
	DDX_Control(pDX, IDC_RDO_XCF, chk_XCF);
	DDX_Control(pDX, IDC_RDO_PCF, chk_PCF);
	DDX_Control(pDX, IDC_PROGRESS, progressBar);
	DDX_Control(pDX, IDC_TXT_BFACT, txt_bfactor);
	DDX_Text(pDX, IDC_TXT_BFACT, s_bfactor);
	DDX_Control(pDX, IDC_TXT_THRESH, txt_thresh);
	DDX_Text(pDX, IDC_TXT_THRESH, s_thresh);
}

BOOL CDMDialog::Create(UINT templateID, CWnd* pParentWnd){return CDialog::Create(IDD, pParentWnd);}

void CDMDialog::OnPaint(){CDialog::OnPaint();}

BOOL CDMDialog::OnEraseBkgnd(CDC* pDC){return CDialog::OnEraseBkgnd(pDC);}

BOOL CDMDialog::OnInitDialog()
{
	CDialog::OnInitDialog();

	int i;

	i = combo_method.AddString("Normal");
	combo_method.SetItemData(i, (DWORD)0);
	i = combo_method.AddString("Overdetermined");
	combo_method.SetItemData(i, (DWORD)1);
	combo_method.SetCurSel(0);
	combo_method.EnableWindow(true);

	chk_XCF.SetCheck(1);
	//chk_PCF.SetCheck(0);

	Devices = OpenCL::GetDeviceList();
	int numDevices = Devices.size();

	if (numDevices < 1)
		{ DMresult << "ERROR: " << "No OpenCL devices found" << DMendl; return false; }

	// Used to index position in list and is stored in combobox
	int ItemData = 0;

	for (std::list<clDevice>::iterator iterator = Devices.begin(), end = Devices.end(); iterator != end; ++iterator)
	{
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

	DMresult << DMendl << "Switching OpenCL device, this can take some time. Please wait..." << DMendl;

	int sel = combo_CLdev.GetCurSel();
	int item = combo_CLdev.GetItemData(sel);

	std::list<clDevice>::iterator iterator = Devices.begin();
	std::advance(iterator, item);
	clDevice dev = *iterator;

	DMresult << "Selected device: " << dev.GetDeviceName() << DMendl;

	clArgStore::SetContext(dev);

	DMresult << "Setting up complete, good to go!" << DMendl;
}

void CDMDialog::OnBnClickedBtnGpa()
{
	if (!(clArgStore::haveDevice))
		{ DMresult << "ERROR: " << " No OpenCL device selected" << DMendl; return; }

	float Bf = boost::lexical_cast<float>(s_bfactor);
	float thresh = boost::lexical_cast<float>(s_thresh);

	int method = combo_method.GetCurSel(); // 0 for normal, 1 for overdetermined
	int corr = chk_PCF.GetCheck(); // 0 for XCF, 1 for PCF

	Align.SetParameters(this, _mtx, Bf, thresh, method, corr);

	Align.StartAlign();
}