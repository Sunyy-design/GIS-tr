
// GIS函数测试程序Dlg.h : 头文件
//

#pragma once


// CGIS函数测试程序Dlg 对话框
class CGIS函数测试程序Dlg : public CDialog
{
// 构造
public:
	CGIS函数测试程序Dlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_GIS_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;
	double d1,f1,m1,x1,d11,f11,m11,x11;
	double d2,f2,m2,x2,d22,f22,m22,x22;
	double Alng;		// A 点经度		角度
	double Alat;		// A 点纬度		角度
	double Blng;		// B 点经度		角度
	double Blat;		// B 点纬度		角度

	long double UTM_AE;
	long double UTM_AN;
	long double UTM_BE;
	long double UTM_BN;

	double dis;			// 距离
	double azi;			// 方位角

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	/*afx_msg void OnBnHotItemChangeRadio1(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnHotItemChangeRadio2(NMHDR *pNMHDR, LRESULT *pResult);*/
	afx_msg void OnBnKillfocusRadio1();
	afx_msg void OnBnSetfocusRadio1();
	afx_msg void OnBnKillfocusRadio2();
	afx_msg void OnBnSetfocusRadio2();
	afx_msg void OnBnClickedRadio1();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButton4();
};
