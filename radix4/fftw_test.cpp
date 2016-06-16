#include <iostream>
#include "fftw3.h"
#include "myfft.h"

using std::cout;
using std::endl;

int main()
{
	size_t batch_size = 1;
	size_t num_size = 8;

	float *x_real, *x_imag;
	float *ref_real, *ref_imag;
	float *test_real, *test_imag;

	x_real = new float[batch_size * num_size];
	x_imag = new float[batch_size * num_size];
	ref_real = new float[batch_size * num_size];
	ref_imag = new float[batch_size * num_size];
	test_real = new float[batch_size * num_size];
	test_imag = new float[batch_size * num_size];

	for (int b = 0; b < batch_size; b++)
	{
		for (int i = 0; i < num_size; i++)
		{
			x_real[b * num_size + i] = static_cast<float>(b * num_size + i);
			x_imag[b * num_size + i] = static_cast<float>(b * num_size + i + 1);
		}
	}

	//fftw
	fftwf_complex *in, *out;
	fftwf_plan plan;
	in = (fftwf_complex*)fftwf_malloc(num_size * sizeof(fftwf_complex));
	out = (fftwf_complex*)fftwf_malloc(num_size * sizeof(fftwf_complex));
	float *in_ptr, *out_ptr;
	in_ptr = (float*)in;
	out_ptr = (float*)out;

	plan = fftwf_plan_dft_1d(num_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	cout << "input vector: " << endl;
	for (int i = 0; i < num_size; i++)
	{
		cout << x_real[i] << ", " << x_imag[i] << " | ";
	}
	cout << endl;
	//execute once per batch
	for (int b = 0; b < batch_size; b++)
	{
		//move the data from x to in; ref to out
		for (int i = 0; i < num_size; i++)
		{
			in_ptr[2 * i] = x_real[b * num_size + i];
			in_ptr[2 * i + 1] = x_imag[b * num_size + i];
		}
		fftwf_execute(plan);
		//move the data from out to ref
		for (int i = 0; i < num_size; i++)
		{
			ref_real[b * num_size + i] = out_ptr[2 * i];
			ref_imag[b * num_size + i] = out_ptr[2 * i + 1];
		}
	}

	cout << "fftw result vector: " << endl;
	for (int i = 0; i < num_size; i++)
	{
		cout << ref_real[i] << ", " << ref_imag[i] << " | ";
	}
	cout << endl;

	fftwf_free(in);
	fftwf_free(out);


	if (num_size == 4)
	{
		for (int b = 0; b < batch_size; b++)
		{
			four_point_fft(x_real + b * 4, x_imag + b * 4, test_real + b * 4, test_imag + b * 4, 0);
		}
	}
	else if (num_size == 8)
	{
		for (int b = 0; b < batch_size; b++)
		{
			eight_point_fft(x_real + b * 4, x_imag + b * 4, test_real + b * 4, test_imag + b * 4);
		}
	}

	cout << "myfft result vector: " << endl;
	for (int i = 0; i < num_size; i++)
	{
		cout << test_real[i] << ", " << test_imag[i] << " | ";
	}
	cout << endl;

	delete[] x_real;
	delete[] x_imag;
	delete[] ref_real;
	delete[] ref_imag;
	delete[] test_real;
	delete[] test_imag;
}