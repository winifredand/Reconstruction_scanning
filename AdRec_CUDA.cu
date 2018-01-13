//The min programming by CUDA7.5;
//By Zengxin
//2018/01/05

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "helper_cuda.h"
#include "helper_functions.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <cstring>
#include <iomanip>

#include <Windows.h>

using namespace std;

const int ReadDataFromFile = 1;
#define Error(s) {printf(s);exit(0);}

__device__ double NewtonF(double a, double b, double c,
	double *x_p, double *y_p, double *z_p, double *q_p)
{
	double R = 1.0;
	double Rd[100], rd[100];
	double cos_theta[100];
	double f_cos[100];
	double ee[100];
	double u[100];
	double L = 0.0;
	double W = 1.0;
	int i;

	for (i = 0; i < 100; i++)
	{
		Rd[i] = (a - x_p[i]) * (a - x_p[i]) + (b - y_p[i]) * (b - y_p[i]) + (c - z_p[i]) * (c - z_p[i]);
		rd[i] = sqrt(Rd[i]);
		cos_theta[i] = ((a - x_p[i]) * x_p[i] + (b - y_p[i]) * y_p[i] + (c - z_p[i]) * z_p[i]) / 2 / rd[i] / R;
		f_cos[i] = 0.999946 + cos_theta[i] * (0.101046 + cos_theta[i] * (-1.040140 + cos_theta[i] * (1.013810 + cos_theta[i] * -0.410953)));
		if (f_cos[i] < 0) f_cos[i] = -f_cos[i];
		ee[i] = exp(-rd[i] / 11);
		u[i] = W * 1 * ee[i] * f_cos[i] / Rd[i];
		L += 2 * ((u[i] - q_p[i]) + log(q_p[i] / u[i]) * q_p[i]);
	}

	return L;
}

void GenerateInitData(double *xx, double *yy, double *zz)
{
	int i;
	double x, y, z;

	srand((unsigned)time(0));

	cout << "The PMT 's position :" << endl;

	for (i = 0; i < 12; i++)
	{
		x = rand() / (double)(RAND_MAX);
		y = rand() / (double)(RAND_MAX);
		z = rand() / (double)(RAND_MAX);

		if ((x * x + y * y + z * z) == 1.0 && x != y || z != y || x != z)
		{
			xx[i] = x;
			yy[i] = y;
			zz[i] = z;
			//printf("%lf	%lf	%lf\n", xx[i], yy[i], zz[i]);
		}
	}

	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 1] = -xx[i];
		yy[i + 12 * 1] = yy[i];
		zz[i + 12 * 1] = zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 2] = xx[i];
		yy[i + 12 * 2] = -yy[i];
		zz[i + 12 * 2] = zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 3] = xx[i];
		yy[i + 12 * 3] = yy[i];
		zz[i + 12 * 3] = -zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 4] = -xx[i];
		yy[i + 12 * 4] = -yy[i];
		zz[i + 12 * 4] = zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 5] = -xx[i];
		yy[i + 12 * 5] = yy[i];
		zz[i + 12 * 5] = -zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 6] = xx[i];
		yy[i + 12 * 6] = -yy[i];
		zz[i + 12 * 6] = -zz[i];
	}
	for (i = 0; i < 12; i++)
	{
		xx[i + 12 * 7] = -xx[i];
		yy[i + 12 * 7] = -yy[i];
		zz[i + 12 * 7] = -zz[i];
	}
	xx[96] = 1.0;
	yy[96] = 0.0;
	zz[96] = 0.0;

	xx[97] = -1.0;
	yy[97] = 0.0;
	zz[97] = 0.0;

	xx[98] = 0.0;
	yy[98] = 1.0;
	zz[98] = 0.0;

	xx[99] = 0.0;
	yy[99] = -1.0;
	zz[99] = 0.0;

}

void GenerateInitDataQ(double *xx, double *yy, double *zz, double *QQ)
{
	int i;
	double Rd[100], rd[100];
	double x = 0.6;
	double y = 0.7;
	double z = 0.2;
	double R = 1.0;
	double cos_theta[100];
	double ee[100];
	double f_cos[100];
	double QQ1[100];
	for (i = 0; i < 100; i++)
	{
		Rd[i] = (x - xx[i]) * (x - xx[i]) + (y - yy[i]) * (y - yy[i]) + (z - zz[i]) * (z - zz[i]);
		rd[i] = sqrt(Rd[i]);
		cos_theta[i] = ((x - xx[i]) * xx[i] + (y - yy[i]) * yy[i] + (z - zz[i]) * zz[i]) / 2 / rd[i] / R;
		f_cos[i] = 0.999946 + cos_theta[i] * (0.101046 + cos_theta[i] * (-1.040140 + cos_theta[i] * (1.013810 + cos_theta[i] * -0.410953)));
		if (f_cos[i] <0) f_cos[i] = -f_cos[i];
		ee[i] = exp(-rd[i] / 11);
		QQ1[i] = ee[i] * f_cos[i] / Rd[i];
		srand((unsigned)time(0));
		QQ[i] = QQ1[i] + rand() / (double)(RAND_MAX);
	}
}

void ReadData(double *x, double *y, double *z, double *Q)   //Read PMT's position and detected charge;
{
	FILE *fp;
	int i;

	if (ReadDataFromFile == 1)
	{
		GenerateInitData(x, y, z);  //Find the features of 100 PMTs that are uniform distributed in ball's surface;

		GenerateInitDataQ(x, y, z, Q);

		for (i = 0; i < 100; i++)
		{
			printf("%lf	%lf	%lf       %lf\n", x[i], y[i], z[i], Q[i]);
		}
		system("pause");
	}
	else
	{
		fp = fopen("data.inp", "r");
		if (fp == NULL) Error("Can not open ctl.inp file. Exit DEM.\n");
		for (i = 0; i < 100; i++)
		{
			fscanf(fp, "%lf	%lf	%lf	%lf", &x[i], &y[i], &z[i], &Q[i]);
		}
		fclose(fp);
	}
}

void CheckData(double *x, double *y, double *z, double *Q)
{
	int i;
	FILE *fp;

	if ((fp = fopen("PMTdata.dat", "w")) == NULL)
		Error("\nCan not open file for debug_particle.\n");

	fprintf(fp, "(x,y,z) E)\n");
	for (i = 0; i < 100; i++)
	{
		fprintf(fp, "(%lf,%lf,%lf)\t", x[i], y[i], z[i]);
		fprintf(fp, "%lf\t", Q[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

__global__ void MINcore(double *d_x_pmt, double *d_y_pmt, double *d_z_pmt, double *d_q_pmt,
	double *d_x, double *d_y, double *d_z, double *d_Like)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	//double like, like_min = 10000000.0;
	//double X[10], Y[10], Z[10];
	//double Like[1000];

//	int l;

	int index, temp, page;
//	int Iter;

	//double down_old = -1, down_new;
	//double delt_old = 0.1, delt_new;

	index = i / 400;
	temp = (i % 400) / 20;
	page = i % 20;

	//double x_new[10], x_old[10];
	//double y_new[10], y_old[10];
	//double z_new[10], z_old[10];
	//double Like_old = 100.0, Like_new;

	//x_old[index] = down_old + delt_old * index;
	//y_old[temp] = down_old + delt_old * temp;
	//z_old[page] = down_old + delt_old * page;

	//for (Iter = 0; Iter < 20; Iter++)
	//{
		//x_new[index] = x_old[index];
		//y_new[temp] = y_old[temp];
		//z_new[page] = z_old[page];

		//Like_new = Like_old;

	d_Like[i] = NewtonF(d_x[index], d_y[temp], d_z[page], d_x_pmt, d_y_pmt, d_z_pmt, d_q_pmt);
	__syncthreads();


		//for (l = 0; l < 1000; l++)
		//{
		//	if (Like_old > Like[l])
		//	{
		//		Like_old = Like[l];
	//	real[0] = x_new[1];
		//real[1] = y_new[10];
	//	real[2] = z_new[1];
		//	}
		//}

		//if (fabs(Like_old - Like_new) < 1E-10) break;
		//if (fabs((Like_old - Like_new) / Like_new) < 1E-10) break;
		/*
		delt_new = delt_old;

		delt_old = delt_new / 5;

		x_old[index] = x_new[index] - delt_new + delt_old * index;
		y_old[temp] = y_new[temp] - delt_new + delt_old * temp;
		z_old[page] = z_new[page] - delt_new + delt_old * page;
	}	*/
}

cudaError_t MINWithCuda(double *x_pmt, double *y_pmt, double *z_pmt, double *q_pmt)
{
	double *d_x_pmt;
	double *d_y_pmt;
	double *d_z_pmt;
	double *d_q_pmt;
	double *d_x;
	double *d_y;
	double *d_z;
	double *d_Like;
	//double *d_real;

	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&d_x_pmt, 100 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&d_y_pmt, 100 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&d_z_pmt, 100 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&d_q_pmt, 100 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&d_x, 20 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&d_y, 20 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&d_z, 20 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&d_Like, 8000 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMemset(d_Like, 0.0, 8000 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}

	cudaStatus = cudaMemcpy(d_x_pmt, x_pmt, 100 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMemcpy(d_y_pmt, y_pmt, 100 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMemcpy(d_z_pmt, z_pmt, 100 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}
	cudaStatus = cudaMemcpy(d_q_pmt, q_pmt, 100 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		goto Error;
	}

	//cudaStatus = cudaMalloc((void**)&d_real, 3 * sizeof(double));
	//if (cudaStatus != cudaSuccess) {
	//	goto Error;
	//}

	double x_new[20], x_old[20];
	double y_new[20], y_old[20];
	double z_new[20], z_old[20];
	double x_down_old = -1, y_down_old = -1, z_down_old = -1;
	double x_down_new, y_down_new, z_down_new;
	double delt_new, delt_old = 0.1;
	double Like_old = 100.0, Like_new;
	double Like_min = 100000000.0;
	double Like[8000];
	int k;
	int IT = 0;
	int temp, index, page;

	for (int i = 0; i < 20; i++)
	{
		x_old[i] = x_down_old + i * delt_old;
		y_old[i] = y_down_old + i * delt_old;
		z_old[i] = z_down_old + i * delt_old;
		printf("%lf	%lf	%lf\n", x_old[i], y_old[i], z_old[i]);
	}

	do
	{
		for (int LL = 0; LL < 20; LL++)
		{
			x_new[LL] = x_old[LL];
			y_new[LL] = y_old[LL];
			z_new[LL] = z_old[LL];
		}

		Like_new = Like_old;

		cudaStatus = cudaMemcpy(d_x, x_new, 20 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_y, y_new, 20 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_z, z_new, 20 * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}

		int Blocks = 8;
		int Threads = 1000;

		MINcore << <Blocks, Threads >> >(d_x_pmt, d_y_pmt, d_z_pmt, d_q_pmt, d_x, d_y, d_z, d_Like);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}

		// Copy output vector from GPU buffer to host memory.
		//cudaStatus = cudaMemcpy(real, d_real, 3 * sizeof(double), cudaMemcpyDeviceToHost);
		//if (cudaStatus != cudaSuccess) {
		//	goto Error;
		//}

		cudaStatus = cudaMemcpy(Like, d_Like, 8000 * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			goto Error;
		}

		for (long i = 0; i < 8000; i++)
		{
			if (Like[i] < Like_min)
			{
				Like_min = Like[i];

				k = i;
			}
		}

		Like_old = Like_min;

		index = k / 400;
		temp = (k % 400) / 20;
		page = k % 20;

		delt_new = delt_old;
		delt_old = delt_new / 10;

		x_down_new = x_down_old;
		y_down_new = y_down_old;
		z_down_new = z_down_old;

		x_down_old = x_new[index] - delt_new;
		y_down_old = y_new[temp] - delt_new;
		z_down_old = z_new[page] - delt_new;

		for (int LL = 0; LL < 20; LL++)
		{
			x_old[LL] = x_down_old + LL * delt_old;
			y_old[LL] = y_down_old + LL * delt_old;
			z_old[LL] = z_down_old + LL * delt_old;
			printf("%.12lf	%.12lf	%.12lf\n", x_old[LL], y_old[LL], z_old[LL]);
		}

		cout << "number of iter = " << IT << endl;
		cout << "Likelihood's value " << setprecision(12) << Like_old << endl;
		cout << "Likelihood's value " << setprecision(12) << Like_new << endl;

		cout << "The x's value " << setprecision(12) << x_new[index] << endl;
		cout << "The y's value " << setprecision(12) << y_new[temp] << endl;
		cout << "The z's value " << setprecision(12) << z_new[page] << endl;

		IT++;

	} while (fabs(Like_old - Like_new) > 1E-20);


Error:
	cudaFree(d_x_pmt);
	cudaFree(d_y_pmt);
	cudaFree(d_z_pmt);
	cudaFree(d_q_pmt);

//	cudaFree(d_real);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_z);
	cudaFree(d_Like);

	return cudaStatus;  
}

int main()
{
	clock_t start, end;  //The time test;
	int i;

	double X_PMT[100];
	double Y_PMT[100];
	double Z_PMT[100];       //The position of PMT;
	double Q_expected[100];  //The expected charge; 

	double X[10];
	double Y[10];
	double Z[10];       //The position of event;

	double x_PMT[2200];
	//double XX[2000];

	ReadData(X_PMT, Y_PMT, Z_PMT, Q_expected);   //Read PMT's position and detected charge;
	CheckData(X_PMT, Y_PMT, Z_PMT, Q_expected);  //Check the random data; 

	for (i = 0; i < 100; i++)
	{
		x_PMT[i * 4] = X_PMT[i];
		x_PMT[i * 4 + 1] = Y_PMT[i];
		x_PMT[i * 4 + 2] = Z_PMT[i];
		x_PMT[i * 4 + 3] = Q_expected[i];
	}
	//double x_real[3];

	start = clock();

	cudaError_t cudaStatus = MINWithCuda(X_PMT, Y_PMT, Z_PMT, Q_expected);
    if (cudaStatus != cudaSuccess) {
        return 1;
    }

	end = clock();

 //   printf("The result is :{%lf,%lf,%lf}\n", x_real[0], x_real[1], x_real[2]);
	cout << "The time is: " << setprecision(12) << (end - start) * CLOCKS_PER_SEC << "s" << endl;

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        return 1;
    }

	system("pause");
    return 0;
}
