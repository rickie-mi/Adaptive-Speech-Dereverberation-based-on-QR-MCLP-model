#ifndef _DSPLIB_H_
#define _DSPLIB_H_
//#include "webrtc/modules/audio_processing/aec/aec_core_internal.h"
#include <math.h>
//DSPF_sp_w_vec说明
//1.x和y数组必须2字节对齐。
//2.nx > 0 且为2的倍数。
//3.y = m * x1 + x2
void DSPF_sp_w_vec(const float *x1, const float *x2, const float m, float *y, const int nx);
//DSPF_sp_vecmul说明
//1.输入输出向量必须2字节对齐。
//2.nx >= 4 且为4的倍数。
void DSPF_sp_vecmul(const float *x1, const float *x2, float *y, const int nx);
void DSPF_sp_mat_mul_cplx(const float *x1, int r1, int c1, const float *x2, int c2, float *y);			//复数矩阵乘法
void DSPF_sp_fftSPxSP(int N, float *ptr_x, float *ptr_w, float *ptr_y, unsigned char *brev, int n_min, int offset, int n_max);
void DSPF_sp_ifftSPxSP(int N, float *ptr_x, float *ptr_w, float *ptr_y, unsigned char *brev, int n_min, int offset, int n_max);
unsigned int bitr(unsigned int a);
void DSPF_sp_vecrecip(const float *x, float *y, const int nx);
float DSPF_sp_vecsum_sq(const float *x, const int nx);
void DSPF_sp_dotp_cplx(const float* x, const float* y, int nx, float *re, float *im);
_Check_return_ __inline float __CRTDECL sqrtsp(_In_ float _X)
{
	return (float)sqrt(_X);
}

_Check_return_ __inline float __CRTDECL cossp(_In_ float _X)
{
	return (float)cos(_X);
}

_Check_return_ __inline float __CRTDECL sinsp(_In_ float _X)
{
	return (float)sin(_X);
}

_Check_return_ __inline float __CRTDECL log10sp(_In_ float _X)
{
	return (float)log10(_X);
}

_Check_return_ __inline float __CRTDECL powsp(_In_ float _X, _In_ float _Y)
{
	return (float)pow(_X, _Y);
}

_Check_return_ __inline float __CRTDECL expsp(_In_ float _X)
{
	return (float)exp(_X);
}

_Check_return_ __inline float __CRTDECL divsp(_In_ float _X, _In_ float _Y)
{
	return (float)(_X / _Y);
}

#endif /* _DSPLIB_H_ */