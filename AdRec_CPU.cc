#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iomanip>

#include <Windows.h>

using namespace std;

#define Error(s) {printf(s);exit(0);}
#define NUM 100
const int ReadDataFromFile = 1;

double Likelihood(double x, double y, double z, double e, double *X, double *Y, double *Z, double *q)
{
	int i;
	double R = 1.0;
	double Rd[100], rd[100];
	double cos_theta[100];
	double f_cos[100];
	double ee[100];
	double u[100];
	double L = 0.0;
	double W = 1.0;

	//cout << "Optical model's parameters" << endl; 
	//printf("i = i\tRd\t        rd\t        cos_theta\t f_cos\t          ee\t        u\t          L\t\n");
	for (i = 0; i < 100; i++)
	{
		Rd[i] = (x - X[i]) * (x - X[i]) + (y - Y[i]) * (y - Y[i]) + (z - Z[i]) * (z - Z[i]);
		rd[i] = sqrt(Rd[i]);
		if (rd[i] <= 0) break;
		cos_theta[i] = ((x - X[i]) * X[i] + (y - Y[i]) * Y[i] + (z - Z[i]) * Z[i]) / 2 / rd[i] / R;
		f_cos[i] = 0.999946 + cos_theta[i] * (0.101046 + cos_theta[i] * (-1.040140 + cos_theta[i] * (1.013810 + cos_theta[i] * -0.410953)));
		if (f_cos[i] <0) f_cos[i] = -f_cos[i];
		ee[i] = exp(-rd[i] / 11);
		u[i] = W * e * ee[i] * f_cos[i] / Rd[i];
		L += 2 * ((u[i] - q[i]) + log(q[i] / u[i]) * q[i]);

		//printf("i = %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", i, Rd[i], rd[i], cos_theta[i], f_cos[i], ee[i], u[i], L);
	}
		//system("pause");
	return L;
}

void GenerateInitData(double *xx, double *yy, double *zz)
{
	int i;
	double x, y, z;

	srand((unsigned)time(0));

	cout << "The PMT 's position :" << endl;

	for (i = 0; i < 100; i++)
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

double MIN(double *x_delt, double *y_delt, double *z_delt, double *xx,
	double *x_pmt, double *y_pmt, double *z_pmt, double *Q_pmt)
{
	double like, like_min = 10000000.0;
	int LL, JJ, KK;
	double x, y, z;
//	int X, Y, Z;
	for (LL = 0; LL < 50; LL++)
	{
		for (KK = 0; KK < 50; KK++)
		{
			for (JJ = 0; JJ < 50; JJ++)
			{
				if (sqrt(x_delt[LL] * x_delt[LL] + y_delt[KK] * y_delt[KK] + z_delt[JJ] * z_delt[JJ]) < 1.0)
				{
					like = Likelihood(x_delt[LL], y_delt[KK], z_delt[JJ], 1.0, x_pmt, y_pmt, z_pmt, Q_pmt);
					//cout << "The " << LL << " " << KK << " " << JJ << "has been done" << endl;
					if (like < like_min)
					{
						like_min = like;
						x = x_delt[LL];
						y = y_delt[KK];
						z = z_delt[JJ];
						//X = LL;
						//Y = KK;
						//Z = LL;
					}
				}

			}
		}
	}

	xx[0] = x;
	xx[1] = y;
	xx[2] = z;

	cout << x << " " << y << " " << z << endl;

	return like_min;
}

int main()
{
	clock_t start, end;
	int k = 0;  //The iteration's frequency;
	double X_PMT[NUM];
	double Y_PMT[NUM];
	double Z_PMT[NUM];       //The position of PMT;
	double Q_expected[NUM];  //The expected charge; 

	ReadData(X_PMT, Y_PMT, Z_PMT, Q_expected);   //Read PMT's position and detected charge;
	CheckData(X_PMT, Y_PMT, Z_PMT, Q_expected);  //Check the random data; 
	

	
	double delt_new, delt_old = 0.04;  //The init step length;
	double x_new[50], x_old[50];
	double y_new[50], y_old[50];
	double z_new[50], z_old[50];
	double Like_new, Like_old = 100.0;
	double x_real[3];
	double x_down_old = -1, y_down_old = -1, z_down_old = -1;
	double x_down_new, y_down_new, z_down_new;
	int LL;

	cout << "x, y, z's init value" << endl;
	for (LL = 0; LL < 50; LL++)
	{
		x_old[LL] = x_down_old + LL * delt_old;
		y_old[LL] = y_down_old + LL * delt_old;
		z_old[LL] = z_down_old + LL * delt_old;
		printf("%lf	%lf	%lf\n", x_old[LL], y_old[LL], z_old[LL]);
	}
	system("pause");

	start = clock();  //The start time;

	do
	{
		for (LL = 0; LL < 50; LL++)
		{
			x_new[LL] = x_old[LL];
			y_new[LL] = y_old[LL];
			z_new[LL] = z_old[LL];
		}
		delt_new = delt_old;
		x_down_new = x_down_old;
		y_down_new = y_down_old;
		z_down_new = z_down_old;

		Like_new = Like_old;
		
		Like_old = MIN(x_new, y_new, z_new, x_real, X_PMT, Y_PMT, Z_PMT, Q_expected);
		cout << x_real[0] << " " << x_real[1] << " " << x_real[2] << endl;

		delt_old = delt_new / 25;
		x_down_old = x_real[0] - delt_new;
		y_down_old = x_real[1] - delt_new;
		z_down_old = x_real[2] - delt_new;

		for (LL = 0; LL < 50; LL++)
		{
			x_old[LL] = x_down_old + LL * delt_old;
			y_old[LL] = y_down_old + LL * delt_old;
			z_old[LL] = z_down_old + LL * delt_old;

			printf("%.12lf	%.12lf	%.12lf\n", x_old[LL], y_old[LL], z_old[LL]);
		}

		cout << "number of iter = " << k << endl;
		cout << "Likelihood's value " << setprecision(12) << Like_old << endl;
		cout << "Likelihood's value " << setprecision(12) << Like_new << endl;

		k++;
	} while (abs(Like_old - Like_new) > 1E-20);  //The eps = 1E-10;

	end = clock();  //The end time;
	cout << "x = " << x_real[0] << endl;
	cout << "y = " << x_real[1] << endl;
	cout << "z = " << x_real[2] << endl;

	cout << "time = " << setprecision(12) << (double)((end - start) / CLOCKS_PER_SEC) << " s" << endl;
	cout << "number of iter = " << k << endl;
	cout << "Likelihood's old value " << setprecision(12) << Like_old << endl;
	cout << "Likelihood's new value " << setprecision(12) << Like_new << endl;

	system("pause");
	return 0;
}
