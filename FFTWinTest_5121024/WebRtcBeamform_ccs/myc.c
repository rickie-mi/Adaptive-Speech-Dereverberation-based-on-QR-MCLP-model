#include "myc.h"

void cnew(mc* mcl, double re, double im)
{
	mcl->re = re;
	mcl->im = im;
}

double cmod(mc a)
{
	double s = sqrt(pow(a.re, 2) + pow(a.im, 2));
	return s;
}


double arctan(double y, double x)
{
	double t;
	if (x > 0)
	{
		if (y > 0)
		{
			t = atan(y / x);
		}
		else
		{
			t = atan(y / x) + 2 * pi;
		}
	}
	else
	{
		t = atan(y / x) + pi;
	}
	return t;
}



mc ccon(mc a)
{
	mc c;
	c.im = -a.im;
	c.re = a.re;
	return c;
}

mc cplus(mc a, mc b)
{
	mc c;
	c.im = a.im + b.im;
	c.re = a.re + b.re;
	return c;
}

mc cminus(mc a, mc b)
{
	mc c;
	c.re = a.re - b.re;
	c.im = a.im - b.im;
	return c;
}

mc cpro(mc a, mc b)
{
	mc c;
	c.re = a.re * b.re - a.im * b.im;
	c.im = a.im * b.re + b.im * a.re;
	return c;
}

mc cdiv(mc a, mc b)
{
	mc bc;
	bc = ccon(b);
	mc c;
	c = cpro(a, bc);
	c.re = c.re / pow(cmod(b), 2);
	c.im = c.im / pow(cmod(b), 2);
	return c;
}

mc cneg(mc a)
{
	mc bc;
	bc.re = -a.re;
	bc.im = -a.im;
	return bc;
}

