
// GIS�������Գ���Dlg.h : ͷ�ļ�
//

#pragma once


// CGIS�������Գ���Dlg �Ի���
class CGIS�������Գ���Dlg : public CDialog
{
// ����
public:
	CGIS�������Գ���Dlg(CWnd* pParent = NULL);	// ��׼���캯��

// �Ի�������
	enum { IDD = IDD_GIS_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��


// ʵ��
protected:
	HICON m_hIcon;
	double d1,f1,m1,x1,d11,f11,m11,x11;
	double d2,f2,m2,x2,d22,f22,m22,x22;
	double Alng;		// A �㾭��		�Ƕ�
	double Alat;		// A ��γ��		�Ƕ�
	double Blng;		// B �㾭��		�Ƕ�
	double Blat;		// B ��γ��		�Ƕ�

	long double UTM_AE;
	long double UTM_AN;
	long double UTM_BE;
	long double UTM_BN;

	double dis;			// ����
	double azi;			// ��λ��

	// ���ɵ���Ϣӳ�亯��
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
