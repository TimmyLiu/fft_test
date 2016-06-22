#ifndef __MY_FFT_H__
#define __MY_FFT_H__

template<typename T>
void four_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag, int stride);

template<typename T>
void eight_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag);

template<typename T>
void three_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag, int stride);

template<typename T>
void nine_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag);
#include "myfft.cpp"

#endif