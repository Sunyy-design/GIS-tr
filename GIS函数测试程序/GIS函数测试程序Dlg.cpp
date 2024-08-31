
// GIS函数测试程序Dlg.cpp : 实现文件
//

#include "stdafx.h"
#include "GIS函数测试程序.h"
#include "GIS函数测试程序Dlg.h"

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define RS_IN					// 输入类参数
#define RS_OUT					// 输出类参数

/****
* 函数功能：将双精度的角度转化为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回弧度
* 2013.7.26 by SaCha
*/
double AngleToRadian( double ang ){
	return M_PI*ang/180.0;
}
/****
* 函数功能：将双精度的弧度转化为角度
* Copyright(C)2013 RabbitSoft Studio
* 返回角度
* 2013.7.26 by SaCha
*/
double RadianToAngle( double rad ){
	return rad*180.0/M_PI;
}

/****
* 函数功能：秒上录入的小数部分转化
* 把一个数123 变成0.123，也就是取出小数点后面的数字
* Copyright(C)2013 RabbitSoft Studio
* 返回小数
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
* 函数功能：录入的四个参数转化为度
* 这个函数有问题 当处理 小数部分为 000123时 前面的0 就都没有了 
* Copyright(C)2013 RabbitSoft Studio
* 返回度数
* 2013.7.26 by SaCha
*/
double FourToAngle( double a ,				//	度
				   double b ,				//	分
				   double c ,				//	秒
				   double d ,				//	秒上的小数部分
				   bool positive)			//	南北true  东西 false
{
	d = dtode(d);
	c += d;
	if( !positive )
		return -1 * ( a + b/60.0 + c/3600.0 );			// 西经 / 南纬
	else
		return a + b/60.0 + c/3600.0;					// 东经 / 北纬
}

/****
* 函数功能：根据经度 求出所在分带的 中央子午线经度
* 单位 度
* Copyright(C)2013 RabbitSoft Studio
* 返回中央子午线经度
* 2013.7.26 by SaCha
********************************
* 为和GE显示一致将 -180划入第1带；180划入第60带；0划入第30带；
*/
double longitude0( double longitude ){

	double longitude0 = 0;

	if( longitude>0 && longitude<180 ){
		// 东经的情况
		longitude0 = ((longitude - fmod(longitude,6.0) + 180)/6 + 1)*6 - 180 -3;
	} else if( longitude>-180 && longitude<0 ){
		// 西经的情况
		longitude0 = ((180 - 6*floor(fabs(longitude)/6) - 6)/6 + 1)*6 - 180 - 3;
	} else if( longitude == 0 ){
		// 0度 第30分带
		longitude0 = -3;
	} else if( longitude == 180 ){
		// 180度 第60分带
		longitude0 = 177;
	} else if( longitude == -180 ){
		// +-180度在一条经线上 第1分带
		longitude0 = -177;
	}

	return longitude0;
}

/****
* 函数功能：计算任意两点距离 参数 x , a 是经度 y ，b是纬度 
* 注意：输入参数为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回两点距离（单位 m）
* R*{arccos[cosb*cosy*cos(a-x)+sinb*siny]}
* 2013.7.26 by SaCha
*/
long double GetDistance(  double lng1 , double lat1 , double lng2 , double lat2  )
{
	// 园的近似算法
	//return R * ( acos( cos(lat2)*cos(lat1)*cos(lng2-lng1) + sin(lat2)*sin(lat1) ) );
	
	// WGS-84 椭球参数
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
* 函数功能：计算任意两点方位角 参数 lng1 , lng2 是经度 lat1 ，lat2是纬度 
* 注意：输入参数为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回两点方位角 （单位 度）
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
* 函数功能：一直起始点坐标、方位角、距离 求终点坐标
* 注意：输入参数为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回未知点经纬度 （单位 ：角度）
* 2013.7.26 by SaCha
*/
void GetTerminus( RS_IN double initlng				// 起始点经度				弧度
						, RS_IN double initlat				// 起始点纬度		弧度
						, RS_IN double azimuth				// 方位角			弧度
						, RS_IN	double dis					// 距离				单位M
						, RS_OUT double &Terlng				// 终点经度			角度
						, RS_OUT double &Terlat )			// 终点纬度			角度
{
	/*
	* 当跨越极点的时候 纬度计算有问题
	*/
	double R = 6378137.0;//6371004.0;

	Terlat = asin( sin(initlat)*cos(dis/R) + cos(initlat)*sin(dis/R)*cos(azimuth) );
	Terlng = initlng + atan2( sin(azimuth)*sin(dis/R)*cos(initlat) , cos(dis/R)-sin(initlat)*sin(Terlat) );

	// 转化成角度
	Terlat = RadianToAngle(Terlat);
	Terlng = RadianToAngle(Terlng);
	return;
}

/****
* 函数功能：已知两点和 方位角 求其交点经纬坐标
* 注意：输入参数为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回未知点经纬度 （单位 ：角度）
* 2013.7.26 by SaCha
*/
void GetIntersection( RS_IN double Alng					// A点经度			弧度
						, RS_IN double Alat				// A点纬度			弧度
						, RS_IN double Aazi				// A方位角			弧度
						, RS_IN double Blng				// B点经度			弧度
						, RS_IN double Blat				// B点纬度			弧度
						, RS_IN double Bazi				// B方位角			弧度
						, RS_OUT double &Terlng			// 交点经度			角度
						, RS_OUT double &Terlat )		// 点纬度			角度
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

	// 转化成角度
	Terlat = RadianToAngle(Terlat);
	Terlng = RadianToAngle(Terlng);
}
/****
* 函数功能：已知四点 求其交点经纬坐标
* 注意：输入参数为弧度
* Copyright(C)2013 RabbitSoft Studio
* 返回未知点经纬度 （单位 ：角度）
* 2013.7.26 by SaCha
*/
void GetIntersection( RS_IN double A1lng					// A1点经度			弧度
						, RS_IN double A1lat				// A1点纬度			弧度
						, RS_IN double A2lng				// A2点经度			弧度
						, RS_IN double A2lat				// A2点纬度			弧度
						, RS_IN double B1lng				// B1点经度			弧度
						, RS_IN double B1lat				// B1点纬度			弧度
						, RS_IN double B2lng				// B2点经度			弧度
						, RS_IN double B2lat				// B2点纬度			弧度
						, RS_OUT double &Terlng1			// 交点1经度			角度
						, RS_OUT double &Terlat1	)		// 交点1纬度			角度
{
	double AziA1_A2 , AziB1_B2;
	// 求出方位角A1->A2 , B1->B2 求1交汇点
	AziA1_A2 = GetAzimuth( A1lng , A1lat , A2lng , A2lat );
	AziB1_B2 = GetAzimuth( B1lng , B1lat , B2lng , B2lat );
	GetIntersection( A1lng , A1lat , AngleToRadian(AziA1_A2) , B1lng , B1lat , AngleToRadian(AziB1_B2) , Terlng1 , Terlat1 );
}

/****
* 函数功能：计算任意经度的子午线弧长
* Copyright(C)2013 RabbitSoft Studio
* 返回子午线弧长结果（单位 m）
* 2013.7.26 by SaCha
* 经过 * 克拉索夫斯基椭球测试数据 符合精度要求
* 克拉索夫斯基于1940年提出的地球椭球，其长半径为6378245米，扁率为1／298.3
* 默认为克氏椭球
*/
long double MeridianArcLength(
							  double lat,				// 大地维度，单位为弧度
							  double a = 6378245,		// 椭球长轴
							  double f =1.0/298.3,		// 椭球扁率
							  int N = 5 )				// 递归项数
{
	int i,n,m;
	long double *k2, e2, ra, c, ff1, k, ff2, sin2;

	k2 = new long double[N];
	for( i=0 ; i<N ; i++ )
		k2[i] = 0.0;

	e2 = f*(2-f);							// 第一偏心率平方
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
* 函数功能：任意椭球的经纬坐标向UTM坐标转换
* Copyright(C)2013 RabbitSoft Studio
* 返回 UTM坐标
* 参数系统默认为 GE的WGS-84椭球参数
* 2013.7.26 by SaCha
*/
void LonAndLatToUTM( RS_OUT long double &UTM_N ,				// 返回的 UTM N数值
					RS_OUT long double &UTM_E ,				// 返回的 UTM E数值
					RS_IN double longitude = 0,			// 已知点经度	角度	单位 度 109.456416 西经为负
					RS_IN double latitude = 0,			// 已知点纬度	角度	单位 度	 35.3516 南纬为负	
					RS_IN double a = 6378137.0,			// 椭球长轴 单位：m
					RS_IN double b = 6356911.94613,		// 椭球短轴 单位：m
					RS_IN double mf = 298.257223563,	// 扁率分母
					RS_IN double e2 = 0.00669437999,		// 第一偏心率平方
					RS_IN double e_2 = 0.0067394967,	// 第二偏心率平方
					RS_IN double k0 = 0.9996	)		// 中央子午线比例系数											
{
	double Rad_lon,Rad_lat;								// 弧度制的经度、纬度
	Rad_lon = AngleToRadian(longitude);
	Rad_lat = AngleToRadian(latitude);

	long double sin1 = AngleToRadian(1.0/3600.0);					// 1秒的 sin值

	// ① 求曲率
	long double f = 1.0/mf;									// 曲率f
	// ② 根据经度求出中央子午线经度，进而取其差值。进而求p
	long double lon0 = longitude0(longitude);				// 所在的中央子午线经度  注意单位是 度
	long double r = (longitude - lon0)*3600;					// 所在经度与中央子午线差值（单位 ：秒）△λ″
	long double p = 0.0001*r;
	// ③ 根据维度求子午线弧长
	long double S = MeridianArcLength( AngleToRadian(latitude) , a , f );
	// ④ 求卯酉线曲率半径
	long double v = a/pow( (1 - e2*sin(Rad_lat)*sin(Rad_lat)) , 0.5 );
	// ⑤ 求T3
	long double T3 = (v*sin(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* (5 - tan(Rad_lat)*tan(Rad_lat) + 9*e_2*cos(Rad_lat)*cos(Rad_lat) + 4*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) ) \
		/ 24.0;
	// ⑥ 求T4
	long double T4 = (v*sin(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* ( 61 - 58*tan(Rad_lat)*tan(Rad_lat) + tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat) + 270*e_2*cos(Rad_lat)*cos(Rad_lat) \
		- 330*tan(Rad_lat)*tan(Rad_lat)*e_2*cos(Rad_lat)*cos(Rad_lat) + 445*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		+ 324*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 680*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		+ 88*e_2*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 600*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 192*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)  ) \
		/ 720.0;
	// ⑦ 求T8
	long double T8 = (v*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0) \
		* ( 5 - 18*tan(Rad_lat)*tan(Rad_lat) + tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat)*tan(Rad_lat) \
		+ 14*e_2*cos(Rad_lat)*cos(Rad_lat) - 58*tan(Rad_lat)*tan(Rad_lat)*e_2*cos(Rad_lat)*cos(Rad_lat) \
		+ 13*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) + 4*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 64*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) \
		- 24*tan(Rad_lat)*tan(Rad_lat)*e_2*e_2*e_2*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat) ) \
		/ 120.0;
	// ⑧ 解N
	UTM_N = S*k0 + (v*sin(Rad_lat)*cos(Rad_lat)*k0)*(sin1*sin1*p*p)*pow(10.0,8.0)/2.0 + T3*(sin1*sin1*sin1*sin1*p*p*p*p)*pow(10.0,16.0) + T4*(sin1*sin1*sin1*sin1*sin1*sin1*p*p*p*p*p*p)*pow(10.0,24.0);
	// ⑨ 解E
	UTM_E = (v*cos(Rad_lat)*k0*sin1*p)*pow(10.0,4.0) + (v*cos(Rad_lat)*cos(Rad_lat)*cos(Rad_lat)*k0)\
		*( 1 - tan(Rad_lat)*tan(Rad_lat) + e_2*cos(Rad_lat)*cos(Rad_lat) )*(sin1*sin1*sin1*p*p*p)*pow(10.0,12.0)/6.0 + T8*(sin1*sin1*sin1*sin1*sin1*p*p*p*p*p)*pow(10.0,20.0) ;

	// ⑩ 为避免负值 加上五百公里
	if( latitude < 0 )
		UTM_N += 10000000;
	UTM_E += 500000;
}


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

	// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	// 实现
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


// CGIS函数测试程序Dlg 对话框




CGIS函数测试程序Dlg::CGIS函数测试程序Dlg(CWnd* pParent /*=NULL*/)
: CDialog(CGIS函数测试程序Dlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CGIS函数测试程序Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CGIS函数测试程序Dlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	/*ON_NOTIFY(BCN_HOTITEMCHANGE, IDC_RADIO1, &CGIS函数测试程序Dlg::OnBnHotItemChangeRadio1)
	ON_NOTIFY(BCN_HOTITEMCHANGE, IDC_RADIO2, &CGIS函数测试程序Dlg::OnBnHotItemChangeRadio2)*/
	ON_BN_KILLFOCUS(IDC_RADIO1, &CGIS函数测试程序Dlg::OnBnKillfocusRadio1)
	ON_BN_SETFOCUS(IDC_RADIO1, &CGIS函数测试程序Dlg::OnBnSetfocusRadio1)
	ON_BN_KILLFOCUS(IDC_RADIO2, &CGIS函数测试程序Dlg::OnBnKillfocusRadio2)
	ON_BN_SETFOCUS(IDC_RADIO2, &CGIS函数测试程序Dlg::OnBnSetfocusRadio2)
	ON_BN_CLICKED(IDC_RADIO1, &CGIS函数测试程序Dlg::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &CGIS函数测试程序Dlg::OnBnClickedRadio2)
	ON_BN_CLICKED(IDC_BUTTON1, &CGIS函数测试程序Dlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CGIS函数测试程序Dlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDOK, &CGIS函数测试程序Dlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON3, &CGIS函数测试程序Dlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CGIS函数测试程序Dlg::OnBnClickedButton4)
END_MESSAGE_MAP()


// CGIS函数测试程序Dlg 消息处理程序

BOOL CGIS函数测试程序Dlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
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

	// 设置此对话框的图标。当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, TRUE);		// 设置小图标

	ShowWindow(SW_SHOWNORMAL);

	// TODO: 在此添加额外的初始化代码

	// 初始化经纬度
	((CComboBox*)GetDlgItem(IDC_COMBO1))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO2))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO3))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO4))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO5))->SetCurSel(0);
	((CComboBox*)GetDlgItem(IDC_COMBO6))->SetCurSel(0);


	// 设置选用度分秒
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(FALSE);//选上
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(TRUE);//不选上
	// 清空控件数据 4 - 11   17 - 24
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

	// 让二组控件只读12 - 15          25 - 28 
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT12))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT14))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT25))->SetReadOnly(0);
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetWindowTextW(TEXT(""));
	((CEdit*)GetDlgItem(IDC_EDIT27))->SetReadOnly(0);


	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CGIS函数测试程序Dlg::OnSysCommand(UINT nID, LPARAM lParam)
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

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CGIS函数测试程序Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CGIS函数测试程序Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


//void CGIS函数测试程序Dlg::OnBnHotItemChangeRadio1(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// 此功能要求 Internet Explorer 6 或更高版本。
//	// 符号 _WIN32_IE 必须是 >= 0x0600。
//	LPNMBCHOTITEM pHotItem = reinterpret_cast<LPNMBCHOTITEM>(pNMHDR);
//	// TODO: 在此添加控件通知处理程序代码
//	*pResult = 0;
//}
//
//void CGIS函数测试程序Dlg::OnBnHotItemChangeRadio2(NMHDR *pNMHDR, LRESULT *pResult)
//{
//	// 此功能要求 Internet Explorer 6 或更高版本。
//	// 符号 _WIN32_IE 必须是 >= 0x0600。
//	LPNMBCHOTITEM pHotItem = reinterpret_cast<LPNMBCHOTITEM>(pNMHDR);
//	// TODO: 在此添加控件通知处理程序代码
//	*pResult = 0;
//}

void CGIS函数测试程序Dlg::OnBnKillfocusRadio1()
{
	// TODO: 在此添加控件通知处理程序代码
}

void CGIS函数测试程序Dlg::OnBnSetfocusRadio1()
{
	// TODO: 在此添加控件通知处理程序代码
}

void CGIS函数测试程序Dlg::OnBnKillfocusRadio2()
{
	// TODO: 在此添加控件通知处理程序代码
}

void CGIS函数测试程序Dlg::OnBnSetfocusRadio2()
{
	// TODO: 在此添加控件通知处理程序代码
}

void CGIS函数测试程序Dlg::OnBnClickedRadio1()
{
	// TODO: 在此添加控件通知处理程序代码
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(TRUE);//选上
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(FALSE);//不选上
	// 清空控件数据 4 - 11   17 - 24
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

	// 让二组控件只读12 - 15          25 - 28 
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

void CGIS函数测试程序Dlg::OnBnClickedRadio2()
{
	// TODO: 在此添加控件通知处理程序代码
	((CButton *)GetDlgItem(IDC_RADIO1))->SetCheck(FALSE);//选上
	((CButton *)GetDlgItem(IDC_RADIO2))->SetCheck(TRUE);//不选上
	// 清空控件数据 4 - 11   17 - 24
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

	// 让二组控件可写12 - 15          25 - 28 
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

void CGIS函数测试程序Dlg::OnBnClickedButton1()
{
	CString str;
	LPCWSTR pstr;	

	// TODO: 在此添加控件通知处理程序代码
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// 录入度分秒形式
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

		// 合并出　Ａ经度　角度
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
		// 合并出　Ａ纬度　角度
		Alat = FourToAngle(d11,f11,m11,x11, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() ) );

	} else {
		// 录入度形式
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT12))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alng = _wtof(pstr);
		str.ReleaseBuffer();
		// 合并出　Ａ经度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() )
			Alng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT14))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alat = _wtof(pstr);
		str.ReleaseBuffer();				
		// 合并出　Ａ纬度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() )
			Alat *= -1;
	}

	// UTM 坐标转换
	::LonAndLatToUTM(UTM_AN,UTM_AE,Alng,Alat);

	str.Format(TEXT("%f"), UTM_AN); 
	((CEdit*)GetDlgItem(IDC_EDIT3))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_AE); 
	((CEdit*)GetDlgItem(IDC_EDIT16))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// 中子午经度显示
	str.Format(TEXT("%f"), ::longitude0(Alng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw1))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS函数测试程序Dlg::OnBnClickedButton2()
{
	// TODO: 在此添加控件通知处理程序代码
	CString str;
	LPCWSTR pstr;

	// TODO: 在此添加控件通知处理程序代码
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// 录入度分秒形式
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

		// 合并出　B经度　角度
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
		// 合并出　B纬度　角度
		Blat = FourToAngle(d22,f22,m22,x22, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() ) );
	} else {
		// 录入度形式
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT25))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blng = _wtof(pstr);
		str.ReleaseBuffer();

		// 合并出　B经度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() )
			Blng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT27))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blat = _wtof(pstr);
		str.ReleaseBuffer();		
		// 合并出　B纬度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() )
			Blat *= -1;
	}

	// UTM 坐标转换
	::LonAndLatToUTM(UTM_BN,UTM_BE,Blng,Blat);

	str.Format(TEXT("%f"), UTM_BN); 
	((CEdit*)GetDlgItem(IDC_EDIT29))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_BE); 
	((CEdit*)GetDlgItem(IDC_EDIT30))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// 中子午经度显示
	str.Format(TEXT("%f"), ::longitude0(Blng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw2))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS函数测试程序Dlg::OnBnClickedOk()
{
	// TODO: 在此添加控件通知处理程序代码
	CString str;
	LPCWSTR pstr;

	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// 录入度分秒形式
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

		// 合并出　Ａ经度　角度
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
		// 合并出　Ａ纬度　角度
		Alat = FourToAngle(d11,f11,m11,x11, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() ) );

		// 录入度分秒形式
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

		// 合并出　B经度　角度
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
		// 合并出　B纬度　角度
		Blat = FourToAngle(d22,f22,m22,x22, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() ) );

	} else {
		// 录入度形式
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT12))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alng = _wtof(pstr);
		str.ReleaseBuffer();
		// 合并出　Ａ经度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO1))->GetCurSel() )
			Alng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT14))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Alat = _wtof(pstr);
		str.ReleaseBuffer();				
		// 合并出　Ａ纬度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO2))->GetCurSel() )
			Alat *= -1;

		// 录入度形式
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT25))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blng = _wtof(pstr);
		str.ReleaseBuffer();

		// 合并出　B经度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO3))->GetCurSel() )
			Blng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT27))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Blat = _wtof(pstr);
		str.ReleaseBuffer();		
		// 合并出　B纬度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO4))->GetCurSel() )
			Blat *= -1;
	}

	// UTM 坐标转换
	::LonAndLatToUTM(UTM_AN,UTM_AE,Alng,Alat);

	str.Format(TEXT("%f"), UTM_AN); 
	((CEdit*)GetDlgItem(IDC_EDIT3))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_AE); 
	((CEdit*)GetDlgItem(IDC_EDIT16))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// UTM 坐标转换
	::LonAndLatToUTM(UTM_BN,UTM_BE,Blng,Blat);

	str.Format(TEXT("%f"), UTM_BN); 
	((CEdit*)GetDlgItem(IDC_EDIT29))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), UTM_BE); 
	((CEdit*)GetDlgItem(IDC_EDIT30))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// 中子午经度显示
	str.Format(TEXT("%f"), ::longitude0(Alng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw1))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// 中子午经度显示
	str.Format(TEXT("%f"), ::longitude0(Blng)); 
	((CEdit*)GetDlgItem(IDC_EDIT_zw2))->SetWindowTextW(str);
	str.ReleaseBuffer();

	// 计算两点距离
	dis = ::GetDistance( ::AngleToRadian(Alng) , ::AngleToRadian(Alat) , ::AngleToRadian(Blng) , ::AngleToRadian(Blat) );
	str.Format(TEXT("%f"), dis); 
	((CEdit*)GetDlgItem(IDC_EDIT1))->SetWindowTextW(str);
	str.ReleaseBuffer();
	// 计算两点方位角
	azi = ::GetAzimuth( ::AngleToRadian(Alng) , ::AngleToRadian(Alat) , ::AngleToRadian(Blng) , ::AngleToRadian(Blat) );
	str.Format(TEXT("%f"), azi); 
	((CEdit*)GetDlgItem(IDC_EDIT2))->SetWindowTextW(str);
	str.ReleaseBuffer();
	//OnOK();
}

void CGIS函数测试程序Dlg::OnBnClickedButton3()
{
	// TODO: 在此添加控件通知处理程序代码
	double Olng,Olat;			// 已知点经纬
	double t_L,t_azi;			// 距离和 方位角
	double d3,f3,m3,x3,d33,f33,m33,x33;
	CString str;
	LPCWSTR pstr;

	// 读取已知点经纬度
	if( 1 == ((CButton *)GetDlgItem(IDC_RADIO1))->GetCheck() )
	{
		// 度分秒
		// 录入度分秒形式
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

		// 合并出　Ａ经度　角度
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
		// 合并出　Ａ纬度　角度
		Olat = FourToAngle(d33,f33,m33,x33, ( 0 == ((CComboBox*)GetDlgItem(IDC_COMBO6))->GetCurSel() ) );
	} else {
		// 度
		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT38))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Olng = _wtof(pstr);
		str.ReleaseBuffer();
		// 合并出　Ａ经度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO5))->GetCurSel() )
			Olng *= -1;

		pstr=TEXT("");
		((CEdit*)GetDlgItem(IDC_EDIT39))->GetWindowTextW(str);		
		pstr = str.GetBuffer();			// CString to double
		Olat = _wtof(pstr);
		str.ReleaseBuffer();				
		// 合并出　Ａ纬度　角度
		if( 1 == ((CComboBox*)GetDlgItem(IDC_COMBO6))->GetCurSel() )
			Olat *= -1;
	}
	// 读距离
	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT40))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	t_L = _wtof(pstr);
	// 读方位角
	pstr=TEXT("");
	((CEdit*)GetDlgItem(IDC_EDIT41))->GetWindowTextW(str);		
	pstr = str.GetBuffer();			// CString to double
	t_azi = _wtof(pstr);

	// 求点
	double TerLng,TerLat;
	::GetTerminus(::AngleToRadian(Olng),::AngleToRadian(Olat),::AngleToRadian(t_azi),t_L,TerLng,TerLat);
	// 显示
	str.Format(TEXT("%f"), TerLng); 
	((CEdit*)GetDlgItem(IDC_EDIT42))->SetWindowTextW(str);
	str.ReleaseBuffer();

	str.Format(TEXT("%f"), TerLat); 
	((CEdit*)GetDlgItem(IDC_EDIT43))->SetWindowTextW(str);
	str.ReleaseBuffer();
}

void CGIS函数测试程序Dlg::OnBnClickedButton4()
{
	// TODO: 在此添加控件通知处理程序代码
	double A1lng,A1lat,A2lng,A2lat;
	double B1lng,B1lat,B2lng,B2lat;
	double Terlng1,Terlat1;

	CString str;
	LPCWSTR pstr;

	// 获取输入参数
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

	// 计算结果而后输出
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
