
// GIS�������Գ���Dlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "GIS�������Գ���.h"
#include "GIS�������Գ���Dlg.h"

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define RS_IN					// ���������
#define RS_OUT					// ��������

/****
* �������ܣ���˫���ȵĽǶ�ת��Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* ���ػ���
* 2013.7.26 by SaCha
*/
double AngleToRadian( double ang ){
	return M_PI*ang/180.0;
}
/****
* �������ܣ���˫���ȵĻ���ת��Ϊ�Ƕ�
* Copyright(C)2013 RabbitSoft Studio
* ���ؽǶ�
* 2013.7.26 by SaCha
*/
double RadianToAngle( double rad ){
	return rad*180.0/M_PI;
}

/****
* �������ܣ�����¼���С������ת��
* ��һ����123 ���0.123��Ҳ����ȡ��С������������
* Copyright(C)2013 RabbitSoft Studio
* ����С��
* 2013.7.26 by SaCha
*/
double dtode( double o )
{
	int z;
	z = (int)o;
	while(1)
	{
		if( z/10 != 0 )
		{
			o = o/10.0;
			z = (int)o;
		}else{
			o = o/10.0;
			break;
		}
	}
	return o;
}

/****
* �������ܣ�¼����ĸ�����ת��Ϊ��
* ������������� ������ С������Ϊ 000123ʱ ǰ���0 �Ͷ�û���� 
* Copyright(C)2013 RabbitSoft Studio
* ���ض���
* 2013.7.26 by SaCha
*/
double FourToAngle( double a ,				//	��
				   double b ,				//	��
				   double c ,				//	��
				   double d ,				//	���ϵ�С������
				   bool positive)			//	�ϱ�true  ���� false
{
	d = dtode(d);
	c += d;
	if( !positive )
		return -1 * ( a + b/60.0 + c/3600.0 );			// ���� / ��γ
	else
		return a + b/60.0 + c/3600.0;					// ���� / ��γ
}

/****
* �������ܣ����ݾ��� ������ڷִ��� ���������߾���
* ��λ ��
* Copyright(C)2013 RabbitSoft Studio
* �������������߾���
* 2013.7.26 by SaCha
********************************
* Ϊ��GE��ʾһ�½� -180�����1����180�����60����0�����30����
*/
double longitude0( double longitude ){

	double longitude0 = 0;

	if( longitude>0 && longitude<180 ){
		// ���������
		longitude0 = ((longitude - fmod(longitude,6.0) + 180)/6 + 1)*6 - 180 -3;
	} else if( longitude>-180 && longitude<0 ){
		// ���������
		longitude0 = ((180 - 6*floor(fabs(longitude)/6) - 6)/6 + 1)*6 - 180 - 3;
	} else if( longitude == 0 ){
		// 0�� ��30�ִ�
		longitude0 = -3;
	} else if( longitude == 180 ){
		// 180�� ��60�ִ�
		longitude0 = 177;
	} else if( longitude == -180 ){
		// +-180����һ�������� ��1�ִ�
		longitude0 = -177;
	}

	return longitude0;
}

/****
* �������ܣ���������������� ���� x , a �Ǿ��� y ��b��γ�� 
* ע�⣺�������Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* ����������루��λ m��
* R*{arccos[cosb*cosy*cos(a-x)+sinb*siny]}
* 2013.7.26 by SaCha
*/
long double GetDistance(  double lng1 , double lat1 , double lng2 , double lat2  )
{
	// ԰�Ľ����㷨
	//return R * ( acos( cos(lat2)*cos(lat1)*cos(lng2-lng1) + sin(lat2)*sin(lat1) ) );
	
	// WGS-84 �������
	double a = 6378137.0;
	double b = 6356752.314245;
	double f = 1/298.257223563;

	double L = lng2 - lng1;
	double U1 = atan( (1-f)*tan(lat1) );
	double U2 = atan( (1-f)*tan(lat2) );

	double sinU1 = sin(U1);
	double sinU2 = sin(U2);
	double cosU1 = cos(U1);
	double cosU2 = cos(U2);

	double lambda = L , lambdaP , iterLimit = 100;

	double sinLambda , cosLambda , sinSigma , cosSigma , sigma , sinAlpha , cosSqAlpha , cos2SigmaM , C , uSq , A , B , deltaSigma , s;
	
	do
	{
		sinLambda = sin(lambda);
		cosLambda = cos(lambda);
		sinSigma = sqrt( (cosU2*sinLambda)*(cosU2*sinLambda) + (cosU1*sinU2-sinU1*cosU2*cosLambda)*(cosU1*sinU2-sinU1*cosU2*cosLambda) );
		
		if( sinSigma == 0 )
			return 0;// co-incident points

		cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
		sigma = atan2( sinSigma , cosSigma );
		sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
		cosSqAlpha = 1 - sinAlpha*sinAlpha;
		cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
		//if( isnan(cos2SigmaM) ) cos2SigmaM = 0;
		C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
		lambdaP = lambda;
		lambda = L + ( 1 - C )*f*sinAlpha*( sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*( -1+2*cos2SigmaM*cos2SigmaM ) ) );
	} while ( abs( lambda - lambdaP ) > 1e-12 && --iterLimit>0 );

	if(iterLimit == 0) 
		return -1;

	uSq = cosSqAlpha*(a*a - b*b)/(b*b);
	A = 1 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)));
	B = uSq/1024.0 * (256.0+uSq*(-128.0+uSq*(74.0-47.0*uSq)));
	deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
	s = b*A*(sigma-deltaSigma);
	return s;
}

/****
* �������ܣ������������㷽λ�� ���� lng1 , lng2 �Ǿ��� lat1 ��lat2��γ�� 
* ע�⣺�������Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* �������㷽λ�� ����λ �ȣ�
* 2013.7.26 by SaCha
*/
double GetAzimuth( double lng1 , double lat1 , double lng2 , double lat2 )
{
	double azi;
	azi = atan2( sin(lng2-lng1)*cos(lat2) , cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lng2-lng1) );
	azi = azi/M_PI*180.0;
	if( azi < 0 )
		azi += 360.0;
	return azi;
}
/****
* �������ܣ�һֱ��ʼ�����ꡢ��λ�ǡ����� ���յ�����
* ע�⣺�������Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* ����δ֪�㾭γ�� ����λ ���Ƕȣ�
* 2013.7.26 by SaCha
*/
void GetTerminus( RS_IN double initlng				// ��ʼ�㾭��				����
						, RS_IN double initlat				// ��ʼ��γ��		����
						, RS_IN double azimuth				// ��λ��			����
						, RS_IN	double dis					// ����				��λM
						, RS_OUT double &Terlng				// �յ㾭��			�Ƕ�
						, RS_OUT double &Terlat )			// �յ�γ��			�Ƕ�
{
	/*
	* ����Խ�����ʱ�� γ�ȼ���������
	*/
	double R = 6378137.0;//6371004.0;

	Terlat = asin( sin(initlat)*cos(dis/R) + cos(initlat)*sin(dis/R)*cos(azimuth) );
	Terlng = initlng + atan2( sin(azimuth)*sin(dis/R)*cos(initlat) , cos(dis/R)-sin(initlat)*sin(Terlat) );

	// ת���ɽǶ�
	Terlat = RadianToAngle(Terlat);
	Terlng = RadianToAngle(Terlng);
	return;
}

/****
* �������ܣ���֪����� ��λ�� ���佻�㾭γ����
* ע�⣺�������Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* ����δ֪�㾭γ�� ����λ ���Ƕȣ�
* 2013.7.26 by SaCha
*/
void GetIntersection( RS_IN double Alng					// A�㾭��			����
						, RS_IN double Alat				// A��γ��			����
						, RS_IN double Aazi				// A��λ��			����
						, RS_IN double Blng				// B�㾭��			����
						, RS_IN double Blat				// B��γ��			����
						, RS_IN double Bazi				// B��λ��			����
						, RS_OUT double &Terlng			// ���㾭��			�Ƕ�
						, RS_OUT double &Terlat )		// ��γ��			�Ƕ�
{
	double C12,C21,a1,a2,a3;
	double brngA,brngB;
	double Drtlat = Blat - Alat;
	double Drtlng = Blng - Alng;

	double d12 = 2*asin( sqrt( sin(Drtlat/2.0)*sin(Drtlat/2.0) + cos(Alat)*cos(Blat)*sin(Drtlng/2.0)*sin(Drtlng/2.0) ) ) ;
	
	if( 0 == d12 )
		return;

	brngA = acos( (sin(Blat) - sin(Alat)*cos(d12))/(sin(d12)*cos(Alat)) );
	brngB = acos( (sin(Alat) - sin(Blat)*cos(d12))/(sin(d12)*cos(Blat)) );

	if( sin(Blng-Alng) > 0 )
	{
		C12 = brngA;
		C21 = 2*M_PI - brngB;
	} else {
		C12 = 2*M_PI - brngA;
		C21 = brngB;
	}

	a1 = fmod( Aazi - C12 + M_PI , 2*M_PI ) - M_PI;
	a2 = fmod( C21 - Bazi + M_PI , 2*M_PI ) - M_PI;
	a3 = acos( -1*cos(a1)*cos(a2) + sin(a1)*sin(a2)*cos(d12) );

	double d13 = atan2( sin(d12)*sin(a1)*sin(a2) , cos(a2) + cos(a1)*cos(a3) );
	Terlat = asin( sin(Alat)*cos(d13) + cos(Alat)*sin(d13)*cos(Aazi) );
	double Drtlng13 = atan2( sin(Aazi)*sin(d13)*cos(Alat) , cos(d13)-sin(Alat)*sin(Terlat) );
	Terlng = fmod( Alng + Drtlng13 + M_PI , 2*M_PI ) - M_PI;

	// ת���ɽǶ�
	Terlat = RadianToAngle(Terlat);
	Terlng = RadianToAngle(Terlng);
}
/****
* �������ܣ���֪�ĵ� ���佻�㾭γ����
* ע�⣺�������Ϊ����
* Copyright(C)2013 RabbitSoft Studio
* ����δ֪�㾭γ�� ����λ ���Ƕȣ�
* 2013.7.26 by SaCha
*/
void GetIntersection( RS_IN double A1lng					// A1�㾭��			����
						, RS_IN double A1lat				// A1��γ��			����
						, RS_IN double A2lng				// A2�㾭��			����
						, RS_IN double A2lat				// A2��γ��			����
						, RS_IN double B1lng				// B1�㾭��			����
						, RS_IN double B1lat				// B1��γ��			����
						, RS_IN double B2lng				// B2�㾭��			����
						, RS_IN double B2lat				// B2��γ��			����
						, RS_OUT double &Terlng1			// ����1����			�Ƕ�
						, RS_OUT double &Terlat1	)		// ����1γ��			�Ƕ�
{
	double AziA1_A2 , AziB1_B2;
	// �����λ��A1->A2 , B1->B2 ��1�����
	AziA1_A2 = GetAzimuth( A1lng , A1lat , A2lng , A2lat );
	AziB1_B2 = GetAzimuth( B1lng , B1lat , B2lng , B2lat );
	GetIntersection( A1lng , A1lat , AngleToRadian(AziA1_A2) , B1lng , B1lat , AngleToRadian(AziB1_B2) , Terlng1 , Terlat1 );
}

/****
* �������ܣ��������⾭�ȵ������߻���
* Copyright(C)2013 RabbitSoft Studio
* ���������߻����������λ m��
* 2013.7.26 by SaCha
* ���� * ��������˹������������� ���Ͼ���Ҫ��
* ��������˹����1940������ĵ��������䳤�뾶Ϊ6378245�ף�����Ϊ1��298.3
* Ĭ��Ϊ��������
*/
long double MeridianArcLength(
							  double lat,				// ���ά�ȣ���λΪ����
							  double a = 6378245,		// ������
							  double f =1.0/298.3,		// �������
							  int N = 5 )				// �ݹ�����
{
	int i,n,m;
	long double *k2, e2, ra, c, ff1, k, ff2, sin2;

	k2 = new long double[N];
	for( i=0 ; i<N ; i++ )
		k2[i] = 0.0;

	e2 = f*(2-f);							// ��һƫ����ƽ��
	ra = a*(1-e2);

	for( c=1.0,n=1 ; n<=N ; n++ )
	{
		c *= (2*n-1.0)*(2*n+1.0)/(4*n*n)*e2;
		for( m=0 ; m<n ; m++ )
			k2[m] += c;
	}

	ff1 = 1.0 + k2[0];
	ff2 = -k2[0];
	sin2 = sin(lat)*sin(lat);

	for( k=1.0,n=1 ; n<N ; n++ )
	{
		k *= 2*n/(2*n+1.0)*sin2;
		ff2 += -k2[n]*k;
	}

	return ra*(lat*ff1 + 0.5*ff2*sin(2.0*lat));
}

/****
* �������ܣ���������ľ�γ������UTM����ת��
* Copyright(C)2013 RabbitSoft Studio
* ���� UTM����
* ����ϵͳĬ��Ϊ GE��WGS-84�������
* 2013.7.26 by SaCha
*/
void LonAndLatToUTM( RS_OUT long double &UTM_N ,				// ���ص� UTM N��ֵ
					RS_OUT long double &UTM_E ,				// ���ص� UTM E��ֵ
					RS_IN double longitude = 0,			// ��֪�㾭��	�Ƕ�	��λ �� 109.456416 ����Ϊ��
					RS_IN double latitude = 0,			// ��֪��γ��	�Ƕ�	��λ ��	 35.3516 ��γΪ��	
					RS_IN double a = 6378137.0,			// ������ ��λ��m
					RS_IN double b = 6356911.94613,		// ������� ��λ��m
					RS_IN double mf = 298.257223563,	// ���ʷ�ĸ
					RS_IN double e2 = 0.00669437999,		// ��һƫ����ƽ��
					RS_IN double e_2 = 0.0067394967,	// �ڶ�ƫ����ƽ��
					RS_IN double k0 = 0.9996	)		// ���������߱���ϵ��											
{
	double Rad_lon,Rad_lat;								// �����Ƶľ��ȡ�γ��
	Rad_lon = AngleToRadian(longitude);
	Rad_lat = AngleToRadian(latitude);

	long double sin1 = AngleToRadian(1.0/3600.0);					// 1��� sinֵ

	// �� ������
	long double f = 1.0/mf;									// ����f
	// �� ���ݾ���������������߾��ȣ�����ȡ���ֵ��������p
	long double lon0 = longitude0(longitude);				// ���ڵ����������߾���  ע�ⵥλ�� ��
	long double r = (longitude - lon0)*3600;					// ���ھ��������������߲�ֵ����λ ���룩���ˡ�
	long double p = 0.0001*r;
	// �� ����ά���������߻���
	long double S = MeridianArcLength( AngleToRadian(latitude) , a , f );
	// �� ��î�������ʰ뾶
	long double v = a/pow( (1 - e2*sin(Rad_lat)*sin(Rad_lat)) , 0.5 );
	// �� ��T3
	long double T3 = (v*sin(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* (5 - tan(Rad_lat)*tan(Rad_lat) + 9*e_2*cos(Rad_lat)*cos(Rad_lat) + 4*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) ) \
		/ 24.0;
	// �� ��T4
	long double T4 = (v*sin(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* ( 61 - 58*tan(Rad_lat)*tan(Rad_lat) + tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat) + 270*e_2*cos(Rad_lat)*cos(Rad_lat) \
		- 330*tan(Rad_lat)*tan(Rad_lat)*e_2*cos(Rad_lat)*cos(Rad_lat) + 445*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		+ 324*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 680*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		+ 88*e_2*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 600*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 192*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)  ) \
		/ 720.0;
	// �� ��T8
	long double T8 = (v*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* ( 5 - 18*tan(Rad_lat)*tan(Rad_lat) + tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat) \
		+ 14*e_2*cos(Rad_lat)*cos(Rad_lat) - 58*tan(Rad_lat)*tan(Rad_lat)*e_2*cos(Rad_lat)*cos(Rad_lat) \
		+ 13*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) + 4*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 64*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 24*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) ) \
		/ 120.0;
	// �� ��N
	UTM_N = S*k0 + (v*sin(Rad_lat)*cos(Rad_lat)*k0)*(sin1*sin1*p*p)*pow(10.0,8.0)/2.0 + T3*(sin1*sin1*sin1*sin1*p*p*p*p)*pow(10.0,16.0) + T4*(sin1*sin1*sin1*sin1*sin1*sin1*p*p*p*p*p*p)*pow(10.0,24.0);
	// �� ��E
	UTM_E = (v*cos(Rad_lat)*k0*sin1*p)*pow(10.0,4.0) + (v*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0)\
		*( 1 - tan(Rad_lat)*tan(Rad_lat) + e_2*cos(Rad_lat)*cos(Rad_lat) )*(sin1*sin1*sin1*p*p*p)*pow(10.0,12.0)/6.0 + T8*(sin1*sin1*sin1*sin1*sin1*p*p*p*p*p)*pow(10.0,20.0) ;

	// �� Ϊ���⸺ֵ ������ٹ���
	if( latitude < 0 )
		UTM_N += 10000000;
	UTM_E += 500000;
}


// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

	// �Ի�������
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	// ʵ��
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CGIS�������Գ���Dlg �Ի���




CGIS�������Գ���Dlg::CGIS�������Գ���Dlg(CWnd* pParent /*=NULL*/)
: CDialog(CGIS�������Գ���Dlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CGIS�������Գ���Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CGIS�������Գ���Dlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	/*ON_NOTIFY(BCN_HOTITEMCHANGE, IDC_RADIO1, &CGIS�������Գ���Dlg::OnBnHotItemChangeRadio1)
	ON_NOTIFY(BCN_HOTITEMCHANGE, IDC_RADIO2, &CGIS�������Գ���Dlg::OnBnHotItemChangeRadio2)*/
	ON_BN_KILLFOCUS(IDC_RADIO1, &CGIS�������Գ���Dlg::OnBnKillfocusRadio1)
	ON_BN_SETFOCUS(IDC_RADIO1, &CGIS�������Գ���Dlg::OnBnSetfocusRadio1)
	ON_BN_KILLFOCUS(IDC_RADIO2, &CGIS�������Գ���Dlg::OnBnKillfocusRadio2)
	ON_BN_SETFOCUS(IDC_RADIO2, &CGIS�������Գ���Dlg::OnBnSetfocusRadio2)
	ON_BN_CLICKED(IDC_RADIO1, &CGIS�������Գ���Dlg::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &CGIS�������Գ���Dlg::OnBnClickedRadio2)
	ON_BN_CLICKED(IDC_BUTTON1, &CGIS�������Գ���Dlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CGIS�������Գ���Dlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDOK, &CGIS�������Գ���Dlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON3, &CGIS�������Գ���Dlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CGIS�������Գ���Dlg::OnBnClickedButton4)
END_MESSAGE_MAP()


// CGIS�������Գ���Dlg ��Ϣ�������

BOOL CGIS�������Գ���Dlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// ���ô˶Ի����ͼ�ꡣ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, TRUE);		// ����Сͼ��

	ShowWindow(SW_SHOWNORMAL);

	// TODO: �ڴ���Ӷ���ĳ�ʼ������

	// ��ʼ����γ��
	((CComboBox*)GetDlgItem(IDC_COMBO1))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO2))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO3))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO4))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO5))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO6))->SetCurSel(0);


	// ����ѡ�öȷ���
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(FALSE);//ѡ��
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(TRUE);//��ѡ��
	// ��տؼ����� 4 - 11   17 - 24
	((CEdit*)GetDlgItem(IDC_EDIT4))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT17))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT4))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetReadOnly(1);

	((CEdit*)GetDlgItem(IDC_EDIT17))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetReadOnly(1);

	((CEdit*)GetDlgItem(IDC_EDIT28))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT28))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetReadOnly(0);

	// �ö���ؼ�ֻ��12 - 15          25 - 28 
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetReadOnly(0);


	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

void CGIS�������Գ���Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CGIS�������Գ���Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CGIS�������Գ���Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


//void CGIS�������Գ���Dlg::OnBnHotItemChangeRadio1(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// �˹���Ҫ�� Internet Explorer 6 ����߰汾��
//	// ���� _WIN32_IE ������ >= 0x0600��
//	LPNMBCHOTITEM pHotItem = reinterpret_cast<LPNMBCHOTITEM>(pNMHDR);
//	// TODO: �ڴ���ӿؼ�֪ͨ����������
//	*pResult = 0;
//}
//
//void CGIS�������Գ���Dlg::OnBnHotItemChangeRadio2(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// �˹���Ҫ�� Internet Explorer 6 ����߰汾��
//	// ���� _WIN32_IE ������ >= 0x0600��
//	LPNMBCHOTITEM pHotItem = reinterpret_cast<LPNMBCHOTITEM>(pNMHDR);
//	// TODO: �ڴ���ӿؼ�֪ͨ����������
//	*pResult = 0;
//}

void CGIS�������Գ���Dlg::OnBnKillfocusRadio1()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
}

void CGIS�������Գ���Dlg::OnBnSetfocusRadio1()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
}

void CGIS�������Գ���Dlg::OnBnKillfocusRadio2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
}

void CGIS�������Գ���Dlg::OnBnSetfocusRadio2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
}

void CGIS�������Գ���Dlg::OnBnClickedRadio1()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(TRUE);//ѡ��
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(FALSE);//��ѡ��
	// ��տؼ����� 4 - 11   17 - 24
	((CEdit*)GetDlgItem(IDC_EDIT4))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT17))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT4))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetReadOnly(0);

	((CEdit*)GetDlgItem(IDC_EDIT17))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetReadOnly(0);

	// �ö���ؼ�ֻ��12 - 15          25 - 28 
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetReadOnly(1);

	((CEdit*)GetDlgItem(IDC_EDIT14))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetReadOnly(1);


	((CEdit*)GetDlgItem(IDC_EDIT25))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetReadOnly(1);

	((CEdit*)GetDlgItem(IDC_EDIT27))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetReadOnly(1);

	((CEdit*)GetDlgItem(IDC_EDIT28))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT28))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetReadOnly(1);

}

void CGIS�������Գ���Dlg::OnBnClickedRadio2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(FALSE);//ѡ��
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(TRUE);//��ѡ��
	// ��տؼ����� 4 - 11   17 - 24
	((CEdit*)GetDlgItem(IDC_EDIT4))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT4))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT5))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT6))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT7))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT8))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT9))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT10))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT11))->SetReadOnly(1);


	((CEdit*)GetDlgItem(IDC_EDIT17))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetWindowTextW(TEXT(""));

	((CEdit*)GetDlgItem(IDC_EDIT17))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT18))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT19))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT20))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT21))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT22))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT23))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT24))->SetReadOnly(1);

	// �ö���ؼ���д12 - 15          25 - 28 
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetReadOnly(0);

	((CEdit*)GetDlgItem(IDC_EDIT14))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetReadOnly(0);


	((CEdit*)GetDlgItem(IDC_EDIT25))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetReadOnly(0);

	((CEdit*)GetDlgItem(IDC_EDIT27))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetReadOnly(0);

	((CEdit*)GetDlgItem(IDC_EDIT28))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT28))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT31))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT32))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT33))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT34))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT35))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT36))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT37))->SetReadOnly(1);
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT38))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT39))->SetReadOnly(0);
}

void CGIS�������Գ���Dlg::OnBnClickedButton1()
{
	CString str;
	LPCWSTR pstr;	

	// TODO: �ڴ���ӿؼ�֪ͨ����������
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// ¼��ȷ�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT4))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT5))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT6))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT7))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x1 = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ����������ȡ��Ƕ�
		Alng = FourToAngle(d1,f1,m1,x1, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() ) );

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT8))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT9))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT10))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT11))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x11 = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ�������γ�ȡ��Ƕ�
		Alat = FourToAngle(d11,f11,m11,x11, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() ) );

	} else {
		// ¼�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT12))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alng = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ����������ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() )
			Alng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT14))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alat = _wtof(pstr);
		str.ReleaseBuffer();				
		// �ϲ�������γ�ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() )
			Alat *= -1;
	}

	// UTM ����ת��
	::LonAndLatToUTM(UTM_AN,UTM_AE,Alng,Alat);

	str.Format(TEXT("%f"), UTM_AN); 
	((CEdit*)GetDlgItem(IDC_EDIT3))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_AE); 
	((CEdit*)GetDlgItem(IDC_EDIT16))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// �����羭����ʾ
	str.Format(TEXT("%f"), ::longitude0(Alng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw1))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS�������Գ���Dlg::OnBnClickedButton2()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	CString str;
	LPCWSTR pstr;

	// TODO: �ڴ���ӿؼ�֪ͨ����������
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// ¼��ȷ�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT17))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT18))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT19))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT20))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x2 = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ�����B���ȡ��Ƕ�
		Blng = FourToAngle(d2,f2,m2,x2, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() ) );

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT21))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT22))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT23))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT24))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x22 = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ�����Bγ�ȡ��Ƕ�
		Blat = FourToAngle(d22,f22,m22,x22, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() ) );
	} else {
		// ¼�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT25))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blng = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ�����B���ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() )
			Blng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT27))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blat = _wtof(pstr);
		str.ReleaseBuffer();		
		// �ϲ�����Bγ�ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() )
			Blat *= -1;
	}

	// UTM ����ת��
	::LonAndLatToUTM(UTM_BN,UTM_BE,Blng,Blat);

	str.Format(TEXT("%f"), UTM_BN); 
	((CEdit*)GetDlgItem(IDC_EDIT29))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_BE); 
	((CEdit*)GetDlgItem(IDC_EDIT30))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// �����羭����ʾ
	str.Format(TEXT("%f"), ::longitude0(Blng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw2))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS�������Գ���Dlg::OnBnClickedOk()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	CString str;
	LPCWSTR pstr;

	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// ¼��ȷ�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT4))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT5))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT6))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m1 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT7))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x1 = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ����������ȡ��Ƕ�
		Alng = FourToAngle(d1,f1,m1,x1, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() ) );

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT8))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT9))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT10))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m11 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT11))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x11 = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ�������γ�ȡ��Ƕ�
		Alat = FourToAngle(d11,f11,m11,x11, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() ) );

		// ¼��ȷ�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT17))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT18))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT19))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m2 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT20))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x2 = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ�����B���ȡ��Ƕ�
		Blng = FourToAngle(d2,f2,m2,x2, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() ) );

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT21))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT22))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT23))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m22 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT24))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x22 = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ�����Bγ�ȡ��Ƕ�
		Blat = FourToAngle(d22,f22,m22,x22, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() ) );

	} else {
		// ¼�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT12))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alng = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ����������ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() )
			Alng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT14))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alat = _wtof(pstr);
		str.ReleaseBuffer();				
		// �ϲ�������γ�ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() )
			Alat *= -1;

		// ¼�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT25))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blng = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ�����B���ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() )
			Blng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT27))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blat = _wtof(pstr);
		str.ReleaseBuffer();		
		// �ϲ�����Bγ�ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() )
			Blat *= -1;
	}

	// UTM ����ת��
	::LonAndLatToUTM(UTM_AN,UTM_AE,Alng,Alat);

	str.Format(TEXT("%f"), UTM_AN); 
	((CEdit*)GetDlgItem(IDC_EDIT3))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_AE); 
	((CEdit*)GetDlgItem(IDC_EDIT16))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// UTM ����ת��
	::LonAndLatToUTM(UTM_BN,UTM_BE,Blng,Blat);

	str.Format(TEXT("%f"), UTM_BN); 
	((CEdit*)GetDlgItem(IDC_EDIT29))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_BE); 
	((CEdit*)GetDlgItem(IDC_EDIT30))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// �����羭����ʾ
	str.Format(TEXT("%f"), ::longitude0(Alng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw1))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// �����羭����ʾ
	str.Format(TEXT("%f"), ::longitude0(Blng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw2))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// �����������
	dis = ::GetDistance( ::AngleToRadian(Alng) , ::AngleToRadian(Alat) , ::AngleToRadian(Blng) , ::AngleToRadian(Blat) );
	str.Format(TEXT("%f"), dis); 
	((CEdit*)GetDlgItem(IDC_EDIT1))->SetWindowTextW(str);
	str.ReleaseBuffer();
	// �������㷽λ��
	azi = ::GetAzimuth( ::AngleToRadian(Alng) , ::AngleToRadian(Alat) , ::AngleToRadian(Blng) , ::AngleToRadian(Blat) );
	str.Format(TEXT("%f"), azi); 
	((CEdit*)GetDlgItem(IDC_EDIT2))->SetWindowTextW(str);
	str.ReleaseBuffer();
	//OnOK();
}

void CGIS�������Գ���Dlg::OnBnClickedButton3()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	double Olng,Olat;			// ��֪�㾭γ
	double t_L,t_azi;			// ����� ��λ��
	double d3,f3,m3,x3,d33,f33,m33,x33;
	CString str;
	LPCWSTR pstr;

	// ��ȡ��֪�㾭γ��
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// �ȷ���
		// ¼��ȷ�����ʽ
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT28))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d3 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT31))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f3 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT32))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m3 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT33))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x3 = _wtof(pstr);
		str.ReleaseBuffer();

		// �ϲ����������ȡ��Ƕ�
		Olng = FourToAngle(d3,f3,m3,x3, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO5))->GetCurSel() ) );

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT34))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		d33 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT35))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		f33 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT36))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		m33 = _wtof(pstr);
		str.ReleaseBuffer();

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT37))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		x33 = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ�������γ�ȡ��Ƕ�
		Olat = FourToAngle(d33,f33,m33,x33, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO6))->GetCurSel() ) );
	} else {
		// ��
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT38))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Olng = _wtof(pstr);
		str.ReleaseBuffer();
		// �ϲ����������ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO5))->GetCurSel() )
			Olng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT39))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Olat = _wtof(pstr);
		str.ReleaseBuffer();				
		// �ϲ�������γ�ȡ��Ƕ�
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO6))->GetCurSel() )
			Olat *= -1;
	}
	// ������
	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT40))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	t_L = _wtof(pstr);
	// ����λ��
	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT41))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	t_azi = _wtof(pstr);

	// ���
	double TerLng,TerLat;
	::GetTerminus(::AngleToRadian(Olng),::AngleToRadian(Olat),::AngleToRadian(t_azi),t_L,TerLng,TerLat);
	// ��ʾ
	str.Format(TEXT("%f"), TerLng); 
	((CEdit*)GetDlgItem(IDC_EDIT42))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), TerLat); 
	((CEdit*)GetDlgItem(IDC_EDIT43))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS�������Գ���Dlg::OnBnClickedButton4()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	double A1lng,A1lat,A2lng,A2lat;
	double B1lng,B1lat,B2lng,B2lat;
	double Terlng1,Terlat1;

	CString str;
	LPCWSTR pstr;

	// ��ȡ�������
	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT44))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	A1lng = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT45))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	A1lat = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT46))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	A2lng = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT47))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	A2lat = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT48))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	B1lng = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT49))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	B1lat = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT50))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	B2lng = _wtof(pstr);
	str.ReleaseBuffer();

	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT51))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	B2lat = _wtof(pstr);
	str.ReleaseBuffer();

	// �������������
	GetIntersection( AngleToRadian(A1lng) , AngleToRadian(A1lat) , AngleToRadian(A2lng) \
		, AngleToRadian(A2lat) , AngleToRadian(B1lng) , AngleToRadian(B1lat) , AngleToRadian(B2lng)\
		, AngleToRadian(B2lat) , Terlng1 , Terlat1 );

	str.Format(TEXT("%f"), Terlng1); 
	((CEdit*)GetDlgItem(IDC_EDIT52))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), Terlat1); 
	((CEdit*)GetDlgItem(IDC_EDIT53))->SetWindowTextW(str);
	str.ReleaseBuffer();
}
