//
// FixPointCS
//
// Copyright(c) Jere Sanisalo, Petri Kero
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include <iostream>

#include <cmath>
#include <random>
#include "FixedUtil.h"
#include "Fixed32.h"
#include "Fixed64.h"

#include "UnitTest.h"

void Test32()
{
	double v = 0.0001f;

	while (v < 150.f)
	{
		Fixed32::FP_INT fv = Fixed32::FromDouble(v);
		Fixed32::FP_INT fv_div = Fixed32::DivPrecise(fv, Fixed32::FromDouble(-2.34));
		Fixed32::FP_INT fv_sqrt = Fixed32::Sqrt(fv);
		Fixed32::FP_INT fv_sin = Fixed32::Sin(fv);
		Fixed32::FP_INT fv_rcp = Fixed32::RcpFast(fv);

		std::cout << v
			<< ": div_by_-2.34: " << Fixed32::ToDouble(fv_div)
			<< ", sqrt: " << Fixed32::ToDouble(fv_sqrt)
			<< ", sin: " << Fixed32::ToDouble(fv_sin)
			<< ", rcp: " << Fixed32::ToDouble(fv_rcp)
			<< std::endl;

		// Next number
		v *= 1.5f;
	}
}

void Test64()
{
	double v = 0.0001f;

	while (v < 150.f)
	{
		Fixed64::FP_LONG fv = Fixed64::FromDouble(v);
		Fixed64::FP_LONG fv_div = Fixed64::DivPrecise(fv, Fixed64::FromDouble(-2.34));
		Fixed64::FP_LONG fv_sqrt = Fixed64::Sqrt(fv);
		Fixed64::FP_LONG fv_sin = Fixed64::Sin(fv);
		Fixed64::FP_LONG fv_rcp = Fixed64::RcpFast(fv);
        Fixed64::FP_LONG fv_div2 = Fixed64::Div(fv, Fixed64::FromDouble(-2.34));

		std::cout << v
			<< ": div_by_-2.34: " << Fixed64::ToDouble(fv_div)
            << ": div2_by_-2.34: " << Fixed64::ToDouble(fv_div2)
			<< ", sqrt: " << Fixed64::ToDouble(fv_sqrt)
			<< ", sin: " << Fixed64::ToDouble(fv_sin)
			<< ", rcp: " << Fixed64::ToDouble(fv_rcp)
			<< std::endl;

		// Next number
		v *= 1.5f;
	}
}
namespace cq {

    //////////////////////////////////////////////////////////////////////////
    // Start of 84 -> 02

    static const double pi = 3.14159265358979324;

    //
    // Krasovsky 1940
    //
    // a = 6378245.0, 1/f = 298.3
    // b = a * (1 - f)
    // ee = (a^2 - b^2) / a^2;
    static const double a = 6378245.0;
    static const double ee = 0.00669342162296594323;

    struct PointD
    {
        double x;
        double y;
    };
    inline PointD PointD_make(double x, double y) { PointD o; o.x = x; o.y = y; return o; }
    static bool outOfChina(double lat, double lon)
    {
        if (lon < 72.004 || lon > 137.8347)
            return true;
        if (lat < 0.8293 || lat > 55.8271)
            return true;
        return false;
    }

    static double transformLat(double x, double y)
    {
        double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
        printf("transformLat :: -100 = %f\n", -100.0);
        printf("transformLat :: 2.0 * x = %f\n", 2.0 * x);
        printf("transformLat :: 3.0 * y = %f\n", 3.0 * y);
        printf("transformLat :: 0.2 * y * y = %f\n", 0.2 * y * y);
        printf("transformLat :: 0.1 * x * y = %f\n", 0.1 * x * y);
        printf("transformLat :: 0.2 * sqrt(abs(x)) = %f\n", 0.2 * sqrt(abs(x)));
        printf("transformLat :: ret1 = %f\n\n", ret);
        ret += (20.0 * sin(6.0 * x * pi) + 20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0;

        printf("transformLat :: 40.0 * sin(6.0 * x * pi) = %f\n", 40.0 * sin(6.0 * x * pi));
        printf("transformLat :: 40.0 * sin(2.0 * x * pi) = %f\n", 40.0 * sin(2.0 * x * pi));
        printf("transformLat :: ret2 = %f\n\n", ret);
        ret += (20.0 * sin(y * pi) + 40.0 * sin(y / 3.0 * pi)) * 2.0 / 3.0;
        printf("transformLat :: 40.0 * sin(y * pi) = %f\n", 40.0 * sin(y * pi));
        printf("transformLat :: 80.0 * sin(y / 3.0 * pi) = %f\n", 80.0 * sin(y / 3.0 * pi));
        printf("transformLat :: ret3 = %f\n\n", ret);
        ret += (160.0 * sin(y / 12.0 * pi) + 320 * sin(y * pi / 30.0)) * 2.0 / 3.0;
        printf("transformLat :: 320.0 * sin(y / 12.0 * pi) = %f\n", 320.0 * sin(y / 12.0 * pi));
        printf("transformLat :: 640 * sin(y * pi / 30.0) = %f\n", 640 * sin(y * pi / 30.0));
        printf("transformLat :: ret4 = %f\n\n", ret);
        return ret;
    }

    static double transformLon(double x, double y)
    {
        double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
        printf("transformLon :: 300.0 = %f\n", 300.0);
        printf("transformLon :: x = %f\n", x);
        printf("transformLon :: 2.0 * y = %f\n", 2.0 * y);
        printf("transformLon :: 0.1 * x * x = %f\n", 0.1 * x * x);
        printf("transformLon :: 0.1 * x * y = %f\n", 0.1 * x * y);
        printf("transformLon :: 0.1 * sqrt(abs(x)) = %f\n", 0.1 * sqrt(abs(x)));
        printf("transformLon :: ret1 = %f\n\n", ret);

        ret += (20.0 * sin(6.0 * x * pi) + 20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0;
        printf("transformLon :: 40.0 * sin(6.0 * x * pi) = %f\n", 40.0 * sin(6.0 * x * pi));
        printf("transformLon :: 40.0 * sin(2.0 * x * pi) = %f\n", 40.0 * sin(2.0 * x * pi));
        printf("transformLon :: ret2 = %f\n\n", ret);
        ret += (20.0 * sin(x * pi) + 40.0 * sin(x / 3.0 * pi)) * 2.0 / 3.0;

        printf("transformLon :: 40.0 * sin(x * pi) = %f\n", 40.0 * sin(x * pi));
        printf("transformLon :: 80.0 * sin(x / 3.0 * pi) = %f\n", 80.0 * sin(x / 3.0 * pi));
        printf("transformLon :: ret3 = %f\n\n", ret);
        ret += (150.0 * sin(x / 12.0 * pi) + 300.0 * sin(x / 30.0 * pi)) * 2.0 / 3.0;

        printf("transformLon :: 300.0 * sin(x / 12.0 * pi) = %f\n", 300.0 * sin(x / 12.0 * pi));
        printf("transformLon :: 600 * sin(x / 30.0 * pi) = %f\n", 600 * sin(x / 30.0 * pi));
        printf("transformLon :: ret4 = %f\n\n", ret);
        return ret;
    }

    //
    // World Geodetic System ==> Mars Geodetic System
    PointD Math_wgsToMarsD(PointD wgsCoordinate)
    {
        double dLat;
        double dLon;
        double radLat;
        double magic;
        double sqrtMagic;

        double wgLat = wgsCoordinate.y;
        double wgLon = wgsCoordinate.x;

        if (outOfChina(wgLat, wgLon))
        {
            return wgsCoordinate;
        }

        dLat = transformLat(wgLon - 105.0, wgLat - 35.0);
        dLon = transformLon(wgLon - 105.0, wgLat - 35.0);
        printf("dLat : %f\n", dLat);
        printf("dLon : %f\n", dLon);
        radLat = wgLat / 180.0 * pi;
        printf("wgLat / 180.0 * pi : %f\n", wgLat / 180.0 * pi);
        magic = sin(radLat);
        printf("sin(radLat) : %f\n", sin(radLat));
        magic = 1 - ee * magic * magic;
        printf("1 - ee * magic * magic : %f\n", 1 - ee * magic * magic);
        sqrtMagic = sqrt(magic);
        printf("sqrt(magic) : %f\n", sqrt(magic));
        printf("dLat * 180.0 : %f\n", dLat * 180.0);
        printf("(a * (1 - ee)) / (magic * sqrtMagic) * pi : %f\n", (a * (1 - ee)) / (magic * sqrtMagic) * pi);
        printf("dLon * 180.0 : %f\n", dLon * 180.0);
        printf("(a / sqrtMagic * cos(radLat) * pi) : %f\n", (a / sqrtMagic * cos(radLat) * pi));
        dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
        dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * pi);
        printf("dLat : %f\n", dLat);
        printf("dLon : %f\n", dLon);

        wgsCoordinate.x = wgLon + dLon;
        wgsCoordinate.y = wgLat + dLat;

        return wgsCoordinate;
    }

    PointD Math_marsToWgsD(PointD marsCoordinate)
    {
        PointD p1, p2;
        p1 = marsCoordinate;
        p2 = Math_wgsToMarsD(p1);
        p1.x = marsCoordinate.x - p2.x + p1.x;
        p1.y = marsCoordinate.y - p2.y + p1.y;
        p2 = Math_wgsToMarsD(p1);
        p1.x = marsCoordinate.x - p2.x + p1.x;
        p1.y = marsCoordinate.y - p2.y + p1.y;
        p2 = Math_wgsToMarsD(p1);
        p1.x = marsCoordinate.x - p2.x + p1.x;
        p1.y = marsCoordinate.y - p2.y + p1.y;
        return p1;
    }
}
namespace fixed64 {

    //////////////////////////////////////////////////////////////////////////
    // Start of 84 -> 02

    static const Fixed64::FP_LONG pi = Fixed64::FromDouble(3.14159265358979324);

    //
    // Krasovsky 1940
    //
    // a = 6378245.0, 1/f = 298.3
    // b = a * (1 - f)
    // ee = (a^2 - b^2) / a^2;
    static const Fixed64::FP_LONG a = Fixed64::FromDouble(6378245.0);
    static const Fixed64::FP_LONG ee = Fixed64::FromDouble(0.00669342162296594323);

    struct PointD
    {
        double x;
        double y;
    };
    inline PointD PointD_make(double x, double y) { PointD o; o.x = x; o.y = y; return o; }
    static bool outOfChina(double lat, double lon)
    {
        if (Fixed64::FromDouble(lon) < Fixed64::FromDouble(72.004) || Fixed64::FromDouble(lon) > Fixed64::FromDouble(137.8347))
            return true;
        if (Fixed64::FromDouble(lat) < Fixed64::FromDouble(0.8293) || Fixed64::FromDouble(lat)> Fixed64::FromDouble(55.8271))
            return true;
        return false;
    }

    static Fixed64::FP_LONG transformLat(Fixed64::FP_LONG x, Fixed64::FP_LONG y)
    {
        //double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
        Fixed64::FP_LONG a1 = Fixed64::FromDouble(-100.0);
        Fixed64::FP_LONG a2 = Fixed64::Mul(Fixed64::FromDouble(2.0), x);
        Fixed64::FP_LONG a3 = Fixed64::Mul(Fixed64::FromDouble(3.0), y);
        Fixed64::FP_LONG a4 = Fixed64::Mul(Fixed64::FromDouble(0.2), Fixed64::Mul(y, y));
        Fixed64::FP_LONG a5 = Fixed64::Mul(Fixed64::FromDouble(0.1), Fixed64::Mul(x, y));
        Fixed64::FP_LONG a6 = Fixed64::Mul(Fixed64::FromDouble(0.2), Fixed64::Sqrt(abs(x)));
        Fixed64::FP_LONG ret = Fixed64::Add(a1, Fixed64::Add(a2, Fixed64::Add(a3, Fixed64::Add(a4, Fixed64::Add(a5, a6)))));
        printf("transformLat :: -100.0 = %f\n", Fixed64::ToDouble(a1));
        printf("transformLat :: 2.0 * x = %f\n", Fixed64::ToDouble(a2));
        printf("transformLat :: 3.0 * y = %f\n", Fixed64::ToDouble(a3));
        printf("transformLat :: 0.2 * y * y = %f\n", Fixed64::ToDouble(a4));
        printf("transformLat :: 0.1 * x * y = %f\n", Fixed64::ToDouble(a5));
        printf("transformLat :: 0.2 * sqrt(abs(x)) = %f\n", Fixed64::ToDouble(a6));
        printf("transformLat :: ret1 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (20.0 * sin(6.0 * x * pi) + 40.0 * sin(2.0 * x * pi)) / 3.0;
        Fixed64::FP_LONG b1 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Mul(Fixed64::FromDouble(6.0), x), pi)));
        Fixed64::FP_LONG b2 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Mul(Fixed64::FromDouble(2.0), x), pi)));
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(b1, b2), Fixed64::FromDouble(3.0)));
        printf("transformLat :: 40.0 * sin(6.0 * x * pi) = %f\n", Fixed64::ToDouble(b1));
        printf("transformLat :: 40.0 * sin(2.0 * x * pi) = %f\n", Fixed64::ToDouble(b2));
        printf("transformLat :: ret2 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (40.0 * sin(y * pi) + 80.0 * sin(y / 3.0 * pi)) / 3.0;
        Fixed64::FP_LONG c1 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(y, pi)));
        Fixed64::FP_LONG c2 = Fixed64::Mul(Fixed64::FromDouble(80.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(y, Fixed64::FromDouble(3.0)), pi)));
          
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(c1, c2), Fixed64::FromDouble(3.0)));
        printf("transformLat :: 40.0 * sin(y * pi) = %f\n", Fixed64::ToDouble(c1));
        printf("transformLat :: 80.0 * sin(y / 3.0 * pi) = %f\n", Fixed64::ToDouble(c2));
        printf("transformLat :: ret3 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (320.0 * sin(y / 12.0 * pi) + 640 * sin(y * pi / 30.0)) / 3.0;
        Fixed64::FP_LONG d1 = Fixed64::Mul(Fixed64::FromDouble(320.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(y, Fixed64::FromDouble(12.0)), pi)));
        Fixed64::FP_LONG d2 = Fixed64::Mul(Fixed64::FromDouble(640.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(y, Fixed64::FromDouble(30.0)), pi)));
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(d1, d2), Fixed64::FromDouble(3.0)));

        printf("transformLat :: 320.0 * sin(y / 12.0 * pi) = %f\n", Fixed64::ToDouble(d1));
        printf("transformLat :: 640 * sin(y * pi / 30.0) = %f\n", Fixed64::ToDouble(d2));

        printf("transformLat :: ret4 = %f\n\n", Fixed64::ToDouble(ret));
        return ret;
    }

    static Fixed64::FP_LONG transformLon(Fixed64::FP_LONG x, Fixed64::FP_LONG y)
    {
        //double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
        Fixed64::FP_LONG a1 = Fixed64::FromDouble(300.0);
        Fixed64::FP_LONG a2 = Fixed64::Mul(Fixed64::FromDouble(1.0), x);
        Fixed64::FP_LONG a3 = Fixed64::Mul(Fixed64::FromDouble(2.0), y);
        Fixed64::FP_LONG a4 = Fixed64::Mul(Fixed64::FromDouble(0.1), Fixed64::Mul(x, x));
        Fixed64::FP_LONG a5 = Fixed64::Mul(Fixed64::FromDouble(0.1), Fixed64::Mul(x, y));
        Fixed64::FP_LONG a6 = Fixed64::Mul(Fixed64::FromDouble(0.1), Fixed64::Sqrt(abs(x)));
        Fixed64::FP_LONG ret = Fixed64::Add(a1, Fixed64::Add(a2, Fixed64::Add(a3, Fixed64::Add(a4, Fixed64::Add(a5, a6)))));
        printf("transformLon :: 300.0 = %f\n", Fixed64::ToDouble(a1));
        printf("transformLon :: x = %f\n", Fixed64::ToDouble(a2));
        printf("transformLon :: 2.0 * y = %f\n", Fixed64::ToDouble(a3));
        printf("transformLon :: 0.1 * x * x = %f\n", Fixed64::ToDouble(a4));
        printf("transformLon :: 0.1 * x * y = %f\n", Fixed64::ToDouble(a5));
        printf("transformLon :: 0.1 * sqrt(abs(x)) = %f\n", Fixed64::ToDouble(a6));
        printf("transformLon :: ret1 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (40.0 * sin(6.0 * x * pi) + 40.0 * sin(2.0 * x * pi)) / 3.0;
        Fixed64::FP_LONG b1 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Mul(Fixed64::FromDouble(6.0), x), pi)));
        Fixed64::FP_LONG b2 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Mul(Fixed64::FromDouble(2.0), x), pi)));
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(b1, b2), Fixed64::FromDouble(3.0)));
        printf("transformLon :: 40.0 * sin(6.0 * x * pi) = %f\n", Fixed64::ToDouble(b1));
        printf("transformLon :: 40.0 * sin(2.0 * x * pi) = %f\n", Fixed64::ToDouble(b2));
        printf("transformLon :: ret2 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (40.0 * sin(x * pi) + 80.0 * sin(x / 3.0 * pi)) / 3.0;
        Fixed64::FP_LONG c1 = Fixed64::Mul(Fixed64::FromDouble(40.0), Fixed64::Sin(Fixed64::Mul(x, pi)));
        Fixed64::FP_LONG c2 = Fixed64::Mul(Fixed64::FromDouble(80.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(x, Fixed64::FromDouble(3.0)), pi)));
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(c1, c2), Fixed64::FromDouble(3.0)));

        printf("transformLon :: 40.0 * sin(x * pi) = %f\n", Fixed64::ToDouble(c1));
        printf("transformLon :: 80.0 * sin(x / 3.0 * pi) = %f\n", Fixed64::ToDouble(c2));
        printf("transformLon :: ret3 = %f\n\n", Fixed64::ToDouble(ret));
        //ret += (300.0 * sin(x / 12.0 * pi) + 600 * sin(x / 30.0 * pi)) / 3.0;
        Fixed64::FP_LONG d1 = Fixed64::Mul(Fixed64::FromDouble(300.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(x, Fixed64::FromDouble(12.0)), pi)));
        Fixed64::FP_LONG d2 = Fixed64::Mul(Fixed64::FromDouble(600.0), Fixed64::Sin(Fixed64::Mul(Fixed64::Div(x, Fixed64::FromDouble(30)), pi)));
        ret = Fixed64::Add(ret, Fixed64::Div(Fixed64::Add(d1, d2), Fixed64::FromDouble(3.0)));

        printf("transformLon :: 300.0 * sin(x / 12.0 * pi) = %f\n", Fixed64::ToDouble(d1));
        printf("transformLon :: 600 * sin(x / 30.0 * pi) = %f\n", Fixed64::ToDouble(d2));
        printf("transformLon :: ret4 = %f\n\n", Fixed64::ToDouble(ret));
        return ret;
    }

    //
    // World Geodetic System ==> Mars Geodetic System
    PointD Math_wgsToMarsD(PointD wgsCoordinate)
    {
        Fixed64::FP_LONG dLat;
        Fixed64::FP_LONG dLon;
        Fixed64::FP_LONG radLat;
        Fixed64::FP_LONG magic;
        Fixed64::FP_LONG sqrtMagic;

        double wgLat = wgsCoordinate.y;
        double wgLon = wgsCoordinate.x;

        if (outOfChina(wgLat, wgLon))
        {
            return wgsCoordinate;
        }

        dLat = transformLat(Fixed64::Sub(Fixed64::FromDouble(wgLon), Fixed64::FromDouble(105.0)), Fixed64::Sub(Fixed64::FromDouble(wgLat), Fixed64::FromDouble(35.0)));
        dLon = transformLon(Fixed64::Sub(Fixed64::FromDouble(wgLon), Fixed64::FromDouble(105.0)), Fixed64::Sub(Fixed64::FromDouble(wgLat), Fixed64::FromDouble(35.0)));

        printf("dLat : %f\n", Fixed64::ToDouble(dLat));
        printf("dLon : %f\n", Fixed64::ToDouble(dLon));
        //radLat = wgLat / 180.0 * pi;
        radLat = Fixed64::Mul(Fixed64::Div(Fixed64::FromDouble(wgLat), Fixed64::FromDouble(180.0)), pi);
        //magic = sin(radLat);
        magic = Fixed64::Sin(radLat);
        //magic = 1 - ee * magic * magic;
        magic = Fixed64::FromDouble(1) - Fixed64::Mul(ee,  Fixed64::Mul(magic, magic));
        //sqrtMagic = sqrt(magic);
        sqrtMagic = Fixed64::Sqrt(magic);
        Fixed64::FP_LONG dLatBak = dLat;
        Fixed64::FP_LONG dLonBak = dLon;
        //dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
        dLat = Fixed64::Div(Fixed64::Mul(dLat, Fixed64::FromDouble(180)), Fixed64::Div(Fixed64::Mul(pi, Fixed64::Mul(a, Fixed64::Sub(Fixed64::FromDouble(1), ee))), Fixed64::Mul(magic, sqrtMagic)));
        //dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * pi);
        dLon = Fixed64::Div(Fixed64::Mul(dLon, Fixed64::FromDouble(180)), Fixed64::Mul(Fixed64::Mul(Fixed64::Div(a, sqrtMagic), Fixed64::Cos(radLat)), pi));

        printf("wgLat / 180.0 * pi : %f\n", Fixed64::ToDouble(Fixed64::Mul(Fixed64::Div(Fixed64::FromDouble(wgLat), Fixed64::FromDouble(180.0)), pi)));
       
        printf("sin(radLat) : %f\n", Fixed64::ToDouble(Fixed64::Sin(radLat)));
        printf("1 - ee * magic * magic : %f\n", Fixed64::ToDouble(Fixed64::FromDouble(1) - Fixed64::Mul(ee, Fixed64::Mul(magic, magic))));
        printf("sqrt(magic) : %f\n", Fixed64::ToDouble(Fixed64::Sqrt(magic)));
        printf("dLat * 180.0 : %f\n", Fixed64::ToDouble(Fixed64::Mul(dLatBak, Fixed64::FromDouble(180))));
        printf("(a * (1 - ee)) / (magic * sqrtMagic) * pi : %f\n", Fixed64::ToDouble(Fixed64::Div(Fixed64::Mul(pi, Fixed64::Mul(a, Fixed64::Sub(Fixed64::FromDouble(1), ee))), Fixed64::Mul(magic, sqrtMagic))));
        printf("dLon * 180.0 : %f\n", Fixed64::ToDouble(Fixed64::Mul(dLonBak, Fixed64::FromDouble(180))));
        printf("(a / sqrtMagic * cos(radLat) * pi) : %f\n", Fixed64::ToDouble(Fixed64::Mul(Fixed64::Mul(Fixed64::Div(a, sqrtMagic), Fixed64::Cos(radLat)), pi)));
        printf("dLat : %f\n", Fixed64::ToDouble(dLat));
        printf("dLon : %f\n", Fixed64::ToDouble(dLon));
        wgsCoordinate.x = Fixed64::ToDouble(Fixed64::FromDouble(wgLon) + dLon);
        wgsCoordinate.y = Fixed64::ToDouble(Fixed64::FromDouble(wgLat) + dLat);

        return wgsCoordinate;
    }

    PointD Math_marsToWgsD(PointD marsCoordinate)
    {
        PointD p1, p2;
        p1 = marsCoordinate;
        p2 = Math_wgsToMarsD(p1);
        p1.x = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.x) - Fixed64::FromDouble(p2.x) + Fixed64::FromDouble(p1.x));
        p1.y = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.y) - Fixed64::FromDouble(p2.y) + Fixed64::FromDouble(p1.y));
        p2 = Math_wgsToMarsD(p1);
        p1.x = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.x) - Fixed64::FromDouble(p2.x) + Fixed64::FromDouble(p1.x));
        p1.y = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.y) - Fixed64::FromDouble(p2.y) + Fixed64::FromDouble(p1.y));
        p2 = Math_wgsToMarsD(p1);
        p1.x = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.x) - Fixed64::FromDouble(p2.x) + Fixed64::FromDouble(p1.x));
        p1.y = Fixed64::ToDouble(Fixed64::FromDouble(marsCoordinate.y) - Fixed64::FromDouble(p2.y) + Fixed64::FromDouble(p1.y));
        return p1;
    }
}

int64_t doubleToFixedPoint2(double value, int fractionalBits) {
    if (value > 2147483647) {
        return INT64_MAX;
    }
    if (value < -2147483648) {
        return INT64_MIN;
    }
    // IEEE 754 格式的存储空间
    uint64_t ieee754;

    // 将 double 的内存复制到 uint64_t
    memcpy(&ieee754, &value, sizeof(value));

    // 提取符号位
    int64_t sign = (ieee754 >> 63) & 0x1; // 符号位
    // 提取指数位并偏移
    int64_t exponent = (ieee754 >> 52) & 0x7FF; // 11位指数
    // 提取尾数（Mantissa）
    int64_t mantissa = ieee754 & 0xFFFFFFFFFFFFF; // 52位尾数

    // 指数调整
    if (exponent == 0) {
        // 非规范化数，处理为 0
        return 0;
    }
    else if (exponent == 0x7FF) {
        // 如果指数位全为1，返回无穷大或NAN
        return (sign) ? INT64_MIN : INT64_MAX; // 根据符号返回最大值或最小值
    }

    // 计算有效指数，将指数减去偏移量 (1023) 并加上 32 （因为要乘以 2^32）
    exponent -= 1023;
    int64_t fixedPointValue = ((1LL << 52) + mantissa);
    int rightShift = (52 - (fractionalBits + exponent));
    if (rightShift >= 0) {
        fixedPointValue = fixedPointValue >> rightShift;
    }
    else {
        fixedPointValue = fixedPointValue << (-rightShift);
    }

    // 如果是负数，调整符号
    if (sign) {
        fixedPointValue = -fixedPointValue;
    }
    return fixedPointValue;
}
double fixedPointToDouble2(int64_t fixedPointValue, int fractionalBits) {
    if (fixedPointValue == 0) {
        return 0;
    }
    // Calculate the scaling factor as 2^fractionalBits
    int64_t scalingFactor = 1LL << fractionalBits; // 2^fractionalBits

    // Determine the sign bit, exponent, and mantissa
    uint64_t sign = (fixedPointValue < 0) ? 1 : 0;
    if (sign) {
        fixedPointValue = -fixedPointValue; // Work with positive value for exponent and mantissa
    }
    int64_t mantissaBits = 1LL << 52;
    int64_t exponent = 0;
    if (fixedPointValue < mantissaBits) {
        while (fixedPointValue < mantissaBits)
        {
            fixedPointValue <<= 1;
            exponent--;
        }
        fixedPointValue = fixedPointValue - mantissaBits;
    }
    else {
        int64_t valueGreaterThanScalingFactor = 0;
        int64_t exp = exponent;
        while (fixedPointValue >= mantissaBits)
        {
            exponent = exp;
            valueGreaterThanScalingFactor = fixedPointValue;
            fixedPointValue >>= 1;
            exp++;
        }
        fixedPointValue = valueGreaterThanScalingFactor - mantissaBits;
    }
    exponent += 1023 + 52 - fractionalBits;

    // Calculate the mantissa
    uint64_t mantissa = fixedPointValue & 0xFFFFFFFFFFFFF; // 52 bits for mantissa

    // Combine into IEEE 754 double format
    uint64_t ieee754 = (sign << 63) | ((exponent & 0x7FF) << 52) | (mantissa & 0xFFFFFFFFFFFFF);

    // Convert uint64_t to double using memcpy
    double result;
    std::memcpy(&result, &ieee754, sizeof(double));
    return result;
}
int64_t FromFloat(float v)
{
    // return (FP_LONG)(v * 4294967296.0f);

    int fractionalBits = 32;
    if (v > 2147483647.0f) {
        return INT64_MAX;
    }
    if (v < -2147483648.0f) {
        return INT64_MIN;
    }

    // IEEE 754 格式的存储空间
    uint32_t ieee754;

    // 将 float 的内存复制到 uint32_t
    memcpy(&ieee754, &v, sizeof(v));

    // 提取符号位
    int32_t sign = (ieee754 >> 31) & 0x1; // 符号位
    // 提取指数位并偏移
    int32_t exponent = (ieee754 >> 23) & 0xFF; // 8位指数
    // 提取尾数（Mantissa）
    int32_t mantissa = ieee754 & 0x7FFFFF; // 23位尾数

    // 指数调整
    if (exponent == 0) {
        // 非规范化数，处理为 0
        return 0;
    }
    else if (exponent == 0xFF) {
        // 如果指数位全为1，返回无穷大或NAN
        return (sign) ? INT64_MIN : INT64_MAX; // 根据符号返回最大值或最小值
    }

    // 计算有效指数，将指数减去偏移量 (127) 并加上 32 （因为要乘以 2^32）
    exponent -= 127;
    int64_t fixedPointValue = ((1LL << 23) + mantissa);

    int leftShift = (exponent + fractionalBits - 23LL);
    if (leftShift >= 0) {
        fixedPointValue = fixedPointValue << leftShift;
    }
    else {
        fixedPointValue = fixedPointValue >> (-leftShift);
    }

    // 如果是负数，调整符号
    if (sign) {
        fixedPointValue = -fixedPointValue;
    }

    return fixedPointValue;
}
float ToFloat(int64_t v)
{
    int fractionalBits = 32;
    if (v == 0) {
        return 0.0f;
    }

    // Calculate the scaling factor as 2^fractionalBits
    int64_t scalingFactor = 1LL << fractionalBits; // 2^fractionalBits

    // Determine the sign bit, exponent, and mantissa
    uint32_t sign = (v < 0) ? 1 : 0;
    if (sign) {
        v = -v; // Work with positive value for exponent and mantissa
    }

    int64_t mantissaBits = 1LL << 23; // 23 bits for mantissa in float
    int32_t exponent = 0;

    if (v < mantissaBits) {
        while (v < mantissaBits) {
            v <<= 1;
            exponent--;
        }
        v = v - mantissaBits;
    }
    else {
        int64_t valueGreaterThanScalingFactor = 0;
        int32_t exp = exponent;
        while (v >= mantissaBits) {
            exponent = exp;
            valueGreaterThanScalingFactor = v;
            v >>= 1;
            exp++;
        }
        v = valueGreaterThanScalingFactor - mantissaBits;
    }

    exponent += 127LL + 23LL - fractionalBits; // Adjust exponent for float

    // Calculate the mantissa
    uint32_t mantissa = v & 0x7FFFFF; // 23 bits for mantissa

    // Combine into IEEE 754 float format
    uint32_t ieee754 = (sign << 31) | ((exponent & 0xFF) << 23) | (mantissa & 0x7FFFFF);

    // Convert uint32_t to float using memcpy
    float result;
    memcpy(&result, &ieee754, sizeof(float));
    return result;
}
union FloatUnion {
    float f;
    uint32_t i;
};

int64_t floatToFixedPoint(float num, int shiftBits) {
    FloatUnion fu;
    fu.f = num;

    // 提取符号位
    int sign = (fu.i >> 31) & 1;
    // 提取指数位
    int exponent = (fu.i >> 23) & 0xFF;
    // 提取尾数位
    uint32_t mantissa = fu.i & 0x7FFFFF;

    int64_t fixedPointValue = 0;

    if (exponent == 0) {
        // 非规范化数或零
        if (mantissa == 0) {
            return 0; // 零
        }
        else {
            // 非规范化数
            // 根据 shiftBits 的值决定左移或右移
            if (shiftBits >= 126) {
                fixedPointValue = static_cast<int64_t>(mantissa) << (shiftBits - 126);
            }
            else {
                fixedPointValue = static_cast<int64_t>(mantissa) >> (126 - shiftBits);
            }
        }
    }
    else if (exponent == 0xFF) {
        // 无穷大或NaN
        if (mantissa == 0) {
            // 无穷大
            return sign ? std::numeric_limits<int64_t>::min() : std::numeric_limits<int64_t>::max();
        }
        else {
            // NaN
            return 0; // 或者可以返回一个特定的错误码
        }
    }
    else {
        // 规范化数
        // 计算实际指数
        int actualExponent = exponent - 127;

        if (shiftBits < actualExponent) {
            // shiftBits 小于实际指数，需要右移尾数位
            int rightShift = actualExponent - shiftBits;
            // 右移尾数位并加上隐含的1
            fixedPointValue = (static_cast<int64_t>(0x800000) | static_cast<int64_t>(mantissa)) >> rightShift;
        }
        else {
            // shiftBits 大于或等于实际指数，直接左移
            fixedPointValue = (static_cast<int64_t>(0x800000) | static_cast<int64_t>(mantissa)) << (shiftBits - actualExponent);
        }
    }

    // 应用符号位
    if (sign) {
        fixedPointValue = -fixedPointValue;
    }

    return fixedPointValue;
}
int main()
{
    /*double gg1 = fixedPointToDouble2(doubleToFixedPoint(123.456789123, 16), 16);
    double gg2 = fixedPointToDouble(doubleToFixedPoint(123.456789123, 16), 16);
    double gg3 = fixedPointToDouble2(doubleToFixedPoint2(123.456789123, 16), 16);
    double gg4 = fixedPointToDouble(doubleToFixedPoint2(123.456789123, 16), 16);
    double gg5 = fixedPointToDouble(doubleToFixedPoint2(123.456789123, 32),32);
    double gg6 = fixedPointToDouble2(doubleToFixedPoint2(123.456789123, 32), 32);
    double gg7 = fixedPointToDouble(Fixed64::FromDouble(123.456789123), 32);
    double gg8 = fixedPointToDouble2(Fixed64::FromDouble(123.456789123), 32);*/
    //std::random_device rd;  // 获取一个随机数种子
    //std::mt19937 gen(rd());  // 使用梅森旋转算法的随机数生成器
    //// 生成均匀分布的随机数
    //std::uniform_real_distribution<double> generator(-18000000.0, 18000000.0);
    //double eps = 1e-9;
    //int64_t dd1 = floatToFixedPoint(123.456789, 24);
    //int64_t dd2 = floatToFixedPoint(0.456789, 24);
    //int64_t dd3 = floatToFixedPoint(-123.456789, 24);
    //int64_t dd4 = floatToFixedPoint(-0.456789, 24);
    //for (size_t i = 0; i < 10000000; i++)
    //{
    //    double x = generator(gen);
    //    float y = x;
    //    float yy = ToFloat(FromFloat(y));

    //    if (std::abs(y - yy) > 1e-1) {
    //        printf("y = %f, yy = %f\n", y, yy);
    //    }
    //    int64_t ax = doubleToFixedPoint2(x, 32);
    //    int64_t bx = Fixed64::FromDouble(x);
    //    double axd = fixedPointToDouble2(ax, 32);
    //    double bxd = Fixed64::ToDouble(bx);
    //    if (std::abs(ax - bx) > Fixed64::FromDouble(eps) || std::abs(x - axd) > eps || std::abs(x - bxd) > eps || std::abs(axd - bxd) > eps) {
    //        printf("ax = %lld, bx = %lld, axd = %.9f, bxd = %.9f, x = %.9f\n", ax, bx, axd, bxd, x);
    //    }
    //}
    //cq::PointD cqMarsPoint = cq::Math_wgsToMarsD(cq::PointD_make(127.995483, 29.688841));
    //std::cout << "cqMarsPoint :" << cqMarsPoint.x << "," << cqMarsPoint.y << std::endl;

    //std::cout << "======================" << std::endl;
    //fixed64::PointD fixed64MarsPoint = fixed64::Math_wgsToMarsD(fixed64::PointD_make(127.995483, 29.688841));
    //std::cout << "fixed64MarsPoint :" << fixed64MarsPoint.x << "," << fixed64MarsPoint.y << std::endl;
	std::cout << "Testing 16.16 fixed point numbers.." << std::endl;
	Test32();

	std::cout << std::endl;
	std::cout << "Testing 32.32 fixed point numbers.." << std::endl;
	Test64();

	std::cout << std::endl;
	std::cout << "Executing all unit tests.." << std::endl;
	UnitTest_TestAll();
	std::cout << "Unit tests finished!" << std::endl;

    return 0;
}
