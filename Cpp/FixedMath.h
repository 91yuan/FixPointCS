#ifndef FIXEDMATH_H
#define FIXEDMATH_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct 
{
    double x;
    double y;
} PointD;

PointD PointD_make(double x, double y);
PointD Math_wgsToMarsD(PointD wgsCoordinate);
PointD Math_marsToWgsD(PointD marsCoordinate);
#ifdef __cplusplus
}
#endif
#endif