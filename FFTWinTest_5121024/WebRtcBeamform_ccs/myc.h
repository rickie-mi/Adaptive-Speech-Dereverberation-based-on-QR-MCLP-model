#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define pi 3.1415926535897932	//Բ����

typedef struct {
	double re;
	double im;
} mc;

/*�¸��� myl�޸�*/
void cnew(mc* mcl, double re, double im);   
/*ģ*/
double cmod(mc a);      


/*arctan(y/x)*/
double arctan(double y, double x);
/*����*/
mc ccon(mc a);
/*����a+b*/
mc cplus(mc a, mc b);
/*����a-c*/
mc cminus(mc a, mc b);
/*����a*b*/
mc cpro(mc a, mc b);
/*����a/b*/
mc cdiv(mc a, mc b);
/*����*/
mc cneg(mc a);
