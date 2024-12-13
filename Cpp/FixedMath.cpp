#include "FixedMath.h"
#include "Fixed64.h"
#include <stdio.h>
#define TRUE 1
#define FALSE 0
#define DEBUG_TRACE 0
using namespace Fixed64;
//////////////////////////////////////////////////////////////////////////
    // Start of 84 -> 02

    //const FP_LONG pi = FromDouble(3.14159265358979324);

    //
    // Krasovsky 1940
    //
    // a = 6378245.0, 1/f = 298.3
    // b = a * (1 - f)
    // ee = (a^2 - b^2) / a^2;
    //const FP_LONG a = FromDouble(6378245.0);
    //const FP_LONG ee = FromDouble(0.00669342162296594323);
    const FP_LONG pi = 13493037704LL;//FromDouble(3.14159265358979324);
    const FP_LONG a = 27394353680875520LL;//FromDouble(6378245.0);
    const FP_LONG ee = 28748026LL;//FromDouble(0.00669342162296594323);
    const FP_LONG P_72_004 = 309254825LL;//FromDouble(72.004)
    const FP_LONG P_137_8347 = 591995528753LL;//FromDouble(137.8347)
    const FP_LONG P_0_8293 = 3561816378LL;//FromDouble(0.8293)
    const FP_LONG P_55_8271 = 239775568730LL;//FromDouble(55.8271)
    const FP_LONG N_100 = -429496729600LL;//FromDouble(-100)
    const FP_LONG P_1 = 4294967296LL;//FromDouble(1.0)
    const FP_LONG P_2 = 8589934592LL;//FromDouble(2.0)
    const FP_LONG P_3 = 12884901888LL;//FromDouble(3.0)
    const FP_LONG P_6 = 25769803776LL;//FromDouble(6.0)
    const FP_LONG P_12 = 51539607552LL;//FromDouble(12.0)
    const FP_LONG P_0_2 = 858993459LL;//FromDouble(0.2)
    const FP_LONG P_0_1 = 429496729LL;//FromDouble(0.1)
    const FP_LONG P_30 = 128849018880LL;//FromDouble(30.0)
    const FP_LONG P_35 = 150323855360LL;//FromDouble(35.0)
    const FP_LONG P_40 = 171798691840LL;//FromDouble(40.0)
    const FP_LONG P_80 = 343597383680LL;//FromDouble(80.0)
    const FP_LONG P_105 = 450971566080LL;//FromDouble(105.0)
    const FP_LONG P_180 = 773094113280LL;//FromDouble(180.0)
    const FP_LONG P_300 = 1288490188800LL;//FromDouble(300.0)
    const FP_LONG P_320 = 1374389534720LL;//FromDouble(320.0)
    const FP_LONG P_600 = 2576980377600LL;//FromDouble(600.0)
    const FP_LONG P_640 = 2748779069440LL;//FromDouble(640.0)
    PointD PointD_make(double x, double y) { PointD o; o.x = x; o.y = y; return o; }
    int8_t outOfChina(FP_LONG lat, FP_LONG lon)
    {
        if (lon < P_72_004 || lon > P_137_8347)
            return TRUE;
        if (lat < P_0_8293 || lat > P_55_8271)
            return TRUE;
        return FALSE;
    }

    FP_LONG transformLat(FP_LONG x, FP_LONG y)
    {
        //double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
        FP_LONG a1 = N_100;
        FP_LONG a2 = Mul(P_2, x);
        FP_LONG a3 = Mul(P_3, y);
        FP_LONG a4 = Mul(P_0_2, Mul(y, y));
        FP_LONG a5 = Mul(P_0_1, Mul(x, y));
        FP_LONG a6 = Mul(P_0_2, Sqrt(x < 0 ? -x : x));
        FP_LONG ret = Add(a1, Add(a2, Add(a3, Add(a4, Add(a5, a6)))));
#if DEBUG_TRACE
        printf("transformLat :: -100.0 = %f\n", ToDouble(a1));
        printf("transformLat :: 2.0 * x = %f\n", ToDouble(a2));
        printf("transformLat :: 3.0 * y = %f\n", ToDouble(a3));
        printf("transformLat :: 0.2 * y * y = %f\n", ToDouble(a4));
        printf("transformLat :: 0.1 * x * y = %f\n", ToDouble(a5));
        printf("transformLat :: 0.2 * sqrt(abs(x)) = %f\n", ToDouble(a6));
        printf("transformLat :: ret1 = %f\n\n", ToDouble(ret));
#endif
        //ret += (20.0 * sin(6.0 * x * pi) + 40.0 * sin(2.0 * x * pi)) / 3.0;
        FP_LONG b1 = Mul(P_40, Sin(Mul(Mul(P_6, x), pi)));
        FP_LONG b2 = Mul(P_40, Sin(Mul(Mul(P_2, x), pi)));
        ret = Add(ret, Div(Add(b1, b2), P_3));
#if DEBUG_TRACE
        printf("transformLat :: 40.0 * sin(6.0 * x * pi) = %f\n", ToDouble(b1));
        printf("transformLat :: 40.0 * sin(2.0 * x * pi) = %f\n", ToDouble(b2));
        printf("transformLat :: ret2 = %f\n\n", ToDouble(ret));
#endif
        //ret += (40.0 * sin(y * pi) + 80.0 * sin(y / 3.0 * pi)) / 3.0;
        FP_LONG c1 = Mul(P_40, Sin(Mul(y, pi)));
        FP_LONG c2 = Mul(P_80, Sin(Mul(Div(y, P_3), pi)));

        ret = Add(ret, Div(Add(c1, c2), P_3));
#if DEBUG_TRACE
        printf("transformLat :: 40.0 * sin(y * pi) = %f\n", ToDouble(c1));
        printf("transformLat :: 80.0 * sin(y / 3.0 * pi) = %f\n", ToDouble(c2));
        printf("transformLat :: ret3 = %f\n\n", ToDouble(ret));
#endif
        //ret += (320.0 * sin(y / 12.0 * pi) + 640 * sin(y * pi / 30.0)) / 3.0;
        FP_LONG d1 = Mul(P_320, Sin(Mul(Div(y, P_12), pi)));
        FP_LONG d2 = Mul(P_640, Sin(Mul(Div(y, P_30), pi)));
        ret = Add(ret, Div(Add(d1, d2), P_3));
#if DEBUG_TRACE
        printf("transformLat :: 320.0 * sin(y / 12.0 * pi) = %f\n", ToDouble(d1));
        printf("transformLat :: 640 * sin(y * pi / 30.0) = %f\n", ToDouble(d2));

        printf("transformLat :: ret4 = %f\n\n", ToDouble(ret));
#endif
        return ret;
    }

    FP_LONG transformLon(FP_LONG x, FP_LONG y)
    {
        //double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
        FP_LONG a1 = P_300;
        FP_LONG a2 = Mul(P_1, x);
        FP_LONG a3 = Mul(P_2, y);
        FP_LONG a4 = Mul(P_0_1, Mul(x, x));
        FP_LONG a5 = Mul(P_0_1, Mul(x, y));
        FP_LONG a6 = Mul(P_0_1, Sqrt(x < 0 ? -x : x));
        FP_LONG ret = Add(a1, Add(a2, Add(a3, Add(a4, Add(a5, a6)))));
#if DEBUG_TRACE
        printf("transformLon :: 300.0 = %f\n", ToDouble(a1));
        printf("transformLon :: x = %f\n", ToDouble(a2));
        printf("transformLon :: 2.0 * y = %f\n", ToDouble(a3));
        printf("transformLon :: 0.1 * x * x = %f\n", ToDouble(a4));
        printf("transformLon :: 0.1 * x * y = %f\n", ToDouble(a5));
        printf("transformLon :: 0.1 * sqrt(abs(x)) = %f\n", ToDouble(a6));
        printf("transformLon :: ret1 = %f\n\n", ToDouble(ret));
#endif
        //ret += (40.0 * sin(6.0 * x * pi) + 40.0 * sin(2.0 * x * pi)) / 3.0;
        FP_LONG b1 = Mul(P_40, Sin(Mul(Mul(P_6, x), pi)));
        FP_LONG b2 = Mul(P_40, Sin(Mul(Mul(P_2, x), pi)));
        ret = Add(ret, Div(Add(b1, b2), P_3));
#if DEBUG_TRACE
        printf("transformLon :: 40.0 * sin(6.0 * x * pi) = %f\n", ToDouble(b1));
        printf("transformLon :: 40.0 * sin(2.0 * x * pi) = %f\n", ToDouble(b2));
        printf("transformLon :: ret2 = %f\n\n", ToDouble(ret));
#endif
        //ret += (40.0 * sin(x * pi) + 80.0 * sin(x / 3.0 * pi)) / 3.0;
        FP_LONG c1 = Mul(P_40, Sin(Mul(x, pi)));
        FP_LONG c2 = Mul(P_80, Sin(Mul(Div(x, P_3), pi)));
        ret = Add(ret, Div(Add(c1, c2), P_3));
#if DEBUG_TRACE
        printf("transformLon :: 40.0 * sin(x * pi) = %f\n", ToDouble(c1));
        printf("transformLon :: 80.0 * sin(x / 3.0 * pi) = %f\n", ToDouble(c2));
        printf("transformLon :: ret3 = %f\n\n", ToDouble(ret));
#endif
        //ret += (300.0 * sin(x / 12.0 * pi) + 600 * sin(x / 30.0 * pi)) / 3.0;
        FP_LONG d1 = Mul(P_300, Sin(Mul(Div(x, P_12), pi)));
        FP_LONG d2 = Mul(P_600, Sin(Mul(Div(x, P_30), pi)));
        ret = Add(ret, Div(Add(d1, d2), P_3));
#if DEBUG_TRACE
        printf("transformLon :: 300.0 * sin(x / 12.0 * pi) = %f\n", ToDouble(d1));
        printf("transformLon :: 600 * sin(x / 30.0 * pi) = %f\n", ToDouble(d2));
        printf("transformLon :: ret4 = %f\n\n", ToDouble(ret));
#endif
        return ret;
    }

    //
    // World Geodetic System ==> Mars Geodetic System
    PointD Math_wgsToMarsD(PointD wgsCoordinate)
    {
        FP_LONG dLat;
        FP_LONG dLon;
        FP_LONG radLat;
        FP_LONG magic;
        FP_LONG sqrtMagic;

        FP_LONG wgLat = FromDouble(wgsCoordinate.y);
        FP_LONG wgLon = FromDouble(wgsCoordinate.x);

        if (outOfChina(wgLat, wgLon))
        {
            return wgsCoordinate;
        }

        dLat = transformLat(Sub(wgLon, P_105), Sub(wgLat, P_35));
        dLon = transformLon(Sub(wgLon, P_105), Sub(wgLat, P_35));
#if DEBUG_TRACE
        printf("dLat : %f\n", ToDouble(dLat));
        printf("dLon : %f\n", ToDouble(dLon));
#endif
        //radLat = wgLat / 180.0 * pi;
        radLat = Mul(Div(wgLat, P_180), pi);
        //magic = sin(radLat);
        magic = Sin(radLat);
        //magic = 1 - ee * magic * magic;
        magic = P_1 - Mul(ee, Mul(magic, magic));
        //sqrtMagic = sqrt(magic);
        if (magic < 0) {
            printf("%.11f,%.11f\n", wgsCoordinate.x, wgsCoordinate.y);
        }
        sqrtMagic = Sqrt(magic);
        FP_LONG dLatBak = dLat;
        FP_LONG dLonBak = dLon;
        //dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
        dLat = Div(Mul(dLat, P_180), Div(Mul(pi, Mul(a, Sub(P_1, ee))), Mul(magic, sqrtMagic)));
        //dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * pi);
        dLon = Div(Mul(dLon, P_180), Mul(Mul(Div(a, sqrtMagic), Cos(radLat)), pi));
#if DEBUG_TRACE
        printf("wgLat / 180.0 * pi : %f\n", ToDouble(Mul(Div(wgLat, P_180), pi)));

        printf("sin(radLat) : %f\n", ToDouble(Sin(radLat)));
        printf("1 - ee * magic * magic : %f\n", ToDouble(P_1 - Mul(ee, Mul(magic, magic))));
        printf("sqrt(magic) : %f\n", ToDouble(Sqrt(magic)));
        printf("dLat * 180.0 : %f\n", ToDouble(Mul(dLatBak, P_180)));
        printf("(a * (1 - ee)) / (magic * sqrtMagic) * pi : %f\n", ToDouble(Div(Mul(pi, Mul(a, Sub(P_1, ee))), Mul(magic, sqrtMagic))));
        printf("dLon * 180.0 : %f\n", ToDouble(Mul(dLonBak, P_180)));
        printf("(a / sqrtMagic * cos(radLat) * pi) : %f\n", ToDouble(Mul(Mul(Div(a, sqrtMagic), Cos(radLat)), pi)));
        printf("dLat : %f\n", ToDouble(dLat));
        printf("dLon : %f\n", ToDouble(dLon));
#endif
        wgsCoordinate.x = ToDouble(wgLon + dLon);
        wgsCoordinate.y = ToDouble(wgLat + dLat);

        return wgsCoordinate;
    }

    PointD Math_marsToWgsD(PointD marsCoordinate)
    {
        PointD p1, p2;
        p1 = marsCoordinate;
        p2 = Math_wgsToMarsD(p1);
        p1.x = ToDouble(FromDouble(marsCoordinate.x) - FromDouble(p2.x) + FromDouble(p1.x));
        p1.y = ToDouble(FromDouble(marsCoordinate.y) - FromDouble(p2.y) + FromDouble(p1.y));
        p2 = Math_wgsToMarsD(p1);
        p1.x = ToDouble(FromDouble(marsCoordinate.x) - FromDouble(p2.x) + FromDouble(p1.x));
        p1.y = ToDouble(FromDouble(marsCoordinate.y) - FromDouble(p2.y) + FromDouble(p1.y));
        p2 = Math_wgsToMarsD(p1);
        p1.x = ToDouble(FromDouble(marsCoordinate.x) - FromDouble(p2.x) + FromDouble(p1.x));
        p1.y = ToDouble(FromDouble(marsCoordinate.y) - FromDouble(p2.y) + FromDouble(p1.y));
        return p1;
    }