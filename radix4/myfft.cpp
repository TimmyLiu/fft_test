#define _USE_MATH_DEFINES
#include <math.h>

template<typename T>
void four_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag, int stride)
{
	/*
	y = Fn * x

	for four points
	Fn =  1  1  1  1
	     [1 -i -1  i]
		  1 -1  1 -1
		  1  i -1 -i
	*/

	out_real[0] = in_real[0] + in_real[1 + stride] + in_real[2 + 2 * stride] + in_real[3 + 3 * stride];
	out_imag[0] = in_imag[0] + in_imag[1 + stride] + in_imag[2 + 2 * stride] + in_imag[3 + 3 * stride];

	out_real[1] = in_real[0] + in_imag[1 + stride] - in_real[2 + 2 * stride] - in_imag[3 + 3 * stride];
	out_imag[1] = in_imag[0] - in_real[1 + stride] - in_imag[2 + 2 * stride] + in_real[3 + 3 * stride];

	out_real[2] = in_real[0] - in_real[1 + stride] + in_real[2 + 2 * stride] - in_real[3 + 3 * stride];
	out_imag[2] = in_imag[0] - in_imag[1 + stride] + in_imag[2 + 2 * stride] - in_imag[3 + 3 * stride];

	out_real[3] = in_real[0] - in_imag[1 + stride] - in_real[2 + 2 * stride] + in_imag[3 + 3 * stride];
	out_imag[3] = in_imag[0] + in_real[1 + stride] - in_imag[2 + 2 * stride] - in_real[3 + 3 * stride];
}

template<typename T>
void eight_point_fft(T *in_real, T *in_imag, T *out_real, T *out_imag)
{
	/*
	recursive method
	*/
	int n = 8;
	int m = n / 2;
	//w = exp(-2*pi*i/n) = cos(-2*pi/n) + i*sin(-2*pi/n)
	double w_real = cos(-2 * M_PI / n);
	double w_imag = sin(-2 * M_PI / n);

	double omega_real_0 = 1;
	double omega_imag_0 = 0;

	double omega_real_1 = w_real;
	double omega_imag_1 = w_imag;

	double omega_real_2 = w_real * w_real - w_imag *w_imag;
	double omega_imag_2 = 2 * w_real * w_imag;

	double omega_real_3 = w_real * w_real * w_real - 3 * w_imag * w_imag * w_real;
	double omega_imag_3 = -1 * w_imag * w_imag * w_imag + 3 * w_real * w_real * w_imag;

	T *odd_real, *odd_imag;
	T *even_real, *even_imag;

	odd_real = new T[m];
	odd_imag = new T[m];
	even_real = new T[m];
	even_imag = new T[m];

	four_point_fft(in_real, in_imag, even_real, even_imag, 1);
	four_point_fft(in_real + 1, in_imag + 1, odd_real, odd_imag, 1);

	//twiddle
	T *twiddle_odd_real, *twiddle_odd_imag;
	twiddle_odd_real = new T[m];
	twiddle_odd_imag = new T[m];
	//odd[0] = omega[0] * odd[0]
	twiddle_odd_real[0] = odd_real[0] * omega_real_0 - odd_imag[0] * omega_imag_0;
	twiddle_odd_imag[0] = odd_real[0] * omega_imag_0 + odd_imag[0] * omega_real_0;

	twiddle_odd_real[1] = odd_real[1] * omega_real_1 - odd_imag[1] * omega_imag_1;
	twiddle_odd_imag[1] = odd_real[1] * omega_imag_1 + odd_imag[1] * omega_real_1;

	twiddle_odd_real[2] = odd_real[2] * omega_real_2 - odd_imag[2] * omega_imag_2;
	twiddle_odd_imag[2] = odd_real[2] * omega_imag_2 + odd_imag[2] * omega_real_2;

	twiddle_odd_real[3] = odd_real[3] * omega_real_3 - odd_imag[3] * omega_imag_3;
	twiddle_odd_imag[3] = odd_real[3] * omega_imag_3 + odd_imag[3] * omega_real_3;

	//combine even and twiddle odd
	/*
	y = [im im;im -im]

	y = [1  0  0  0  1  0  0  0]
	    [0  1  0  0  0  1  0  0]
		[0  0  1  0  0  0  1  0]
		[0  0  0  1  0  0  0  1]
		[1  0  0  0 -1  0  0  0]
		[0  1  0  0  0 -1  0  0]
		[0  0  1  0  0  0 -1  0]
		[0  0  0  1  0  0  0 -1]
	
	x = [even_real[0], even_imag[0] ]
	    [even_real[1], even_imag[1] ]
		[even_real[2], even_imag[2] ]
		[even_real[3], even_imag[3] ]
		[twiddle_odd_real[0], twiddle_odd_imag[0] ]
		[twiddle_odd_real[1], twiddle_odd_imag[1] ]
		[twiddle_odd_real[2], twiddle_odd_imag[2] ]
		[twiddle_odd_real[3], twiddle_odd_imag[3] ]
	*/

	out_real[0] = even_real[0] + twiddle_odd_real[0];
	out_imag[0] = even_imag[0] + twiddle_odd_imag[0];

	out_real[1] = even_real[1] + twiddle_odd_real[1];
	out_imag[1] = even_imag[1] + twiddle_odd_imag[1];

	out_real[2] = even_real[2] + twiddle_odd_real[2];
	out_imag[2] = even_imag[2] + twiddle_odd_imag[2];
	
	out_real[3] = even_real[3] + twiddle_odd_real[3];
	out_imag[3] = even_imag[3] + twiddle_odd_imag[3];

	out_real[4] = even_real[0] - twiddle_odd_real[0];
	out_imag[4] = even_imag[0] - twiddle_odd_imag[0];

	out_real[5] = even_real[1] - twiddle_odd_real[1];
	out_imag[5] = even_imag[1] - twiddle_odd_imag[1];

	out_real[6] = even_real[2] - twiddle_odd_real[2];
	out_imag[6] = even_imag[2] - twiddle_odd_imag[2];

	out_real[7] = even_real[3] - twiddle_odd_real[3];
	out_imag[7] = even_imag[3] - twiddle_odd_imag[3];


}