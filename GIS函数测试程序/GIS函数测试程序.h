
// GIS�������Գ���.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CGIS�������Գ���App:
// �йش����ʵ�֣������ GIS�������Գ���.cpp
//

class CGIS�������Գ���App : public CWinAppEx
{
public:
	CGIS�������Գ���App();

// ��д
	public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CGIS�������Գ���App theApp;