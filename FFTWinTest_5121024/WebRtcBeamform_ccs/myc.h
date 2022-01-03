#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define pi 3.1415926535897932	//圆周率

typedef struct {
	double re;
	double im;
} mc;

/*新复数 myl修改*/
void cnew(mc* mcl, double re, double im);   
/*模*/
double cmod(mc a);      


/*arctan(y/x)*/
double arctan(double y, double x);
/*共轭*/
mc ccon(mc a);
/*复数a+b*/
mc cplus(mc a, mc b);
/*复数a-c*/
mc cminus(mc a, mc b);
/*复数a*b*/
mc cpro(mc a, mc b);
/*复数a/b*/
mc cdiv(mc a, mc b);
/*求负数*/
mc cneg(mc a);
