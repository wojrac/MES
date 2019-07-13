#pragma once
double policzWyznacznik(double przedJakobian[2][2])
{
	double wyznacznik = 0;
	wyznacznik = przedJakobian[0][0] * przedJakobian[1][1] - przedJakobian[0][1] * przedJakobian[1][0];
	return wyznacznik;
}
double** odwrocPrzedJakobian(double przedJakobian[2][2])
{
	double** Jakobian = new double*[2];
	for (int i = 0; i < 2; i++)
	{
		Jakobian[i] = new double[2];
	}
	Jakobian[0][0] = przedJakobian[1][1];
	Jakobian[0][1] = (-1)*przedJakobian[0][1];
	Jakobian[1][0] = (-1)*przedJakobian[1][0];
	Jakobian[1][1] = przedJakobian[0][0];
	return Jakobian;
}
double** transponujMnoz(double wektor[4])
{
	double** transpon = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		transpon[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			transpon[i][j] = wektor[i] * wektor[j];
		}
	}
	return transpon;
}
double** liczHpow1(double funkcjeKsztaltu_dS_fun[4][4], double alfaFun, double bokFun)
{
	double** matrixResult = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrixResult[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixResult[i][j] = 0;
		}
	}
	//[pc][fk]
	matrixResult[0][0] = ((funkcjeKsztaltu_dS_fun[0][0] * funkcjeKsztaltu_dS_fun[0][0]) + (funkcjeKsztaltu_dS_fun[1][0] * funkcjeKsztaltu_dS_fun[1][0]))*alfaFun *bokFun / 2;
	matrixResult[0][1] = ((funkcjeKsztaltu_dS_fun[0][0] * funkcjeKsztaltu_dS_fun[0][1]) + (funkcjeKsztaltu_dS_fun[1][0] * funkcjeKsztaltu_dS_fun[1][1]))*alfaFun *bokFun / 2;
	matrixResult[1][0] = ((funkcjeKsztaltu_dS_fun[0][1] * funkcjeKsztaltu_dS_fun[0][0]) + (funkcjeKsztaltu_dS_fun[1][1] * funkcjeKsztaltu_dS_fun[1][0]))*alfaFun *bokFun / 2;
	matrixResult[1][1] = ((funkcjeKsztaltu_dS_fun[0][1] * funkcjeKsztaltu_dS_fun[0][1]) + (funkcjeKsztaltu_dS_fun[1][1] * funkcjeKsztaltu_dS_fun[1][1]))*alfaFun *bokFun / 2;
	return matrixResult;
}
double** liczHpow2(double funkcjeKsztaltu_dS_fun[4][4], double alfaFun, double bokFun)
{
	double** matrixResult = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrixResult[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixResult[i][j] = 0;
		}
	}
	matrixResult[1][1] = ((funkcjeKsztaltu_dS_fun[2][1] * funkcjeKsztaltu_dS_fun[2][1]) + (funkcjeKsztaltu_dS_fun[3][1] * funkcjeKsztaltu_dS_fun[3][1]))*alfaFun *bokFun / 2;
	matrixResult[1][2] = ((funkcjeKsztaltu_dS_fun[2][1] * funkcjeKsztaltu_dS_fun[2][2]) + (funkcjeKsztaltu_dS_fun[3][1] * funkcjeKsztaltu_dS_fun[3][2]))*alfaFun *bokFun / 2;
	matrixResult[2][1] = ((funkcjeKsztaltu_dS_fun[2][2] * funkcjeKsztaltu_dS_fun[2][1]) + (funkcjeKsztaltu_dS_fun[3][2] * funkcjeKsztaltu_dS_fun[3][1]))*alfaFun *bokFun / 2;
	matrixResult[2][2] = ((funkcjeKsztaltu_dS_fun[2][2] * funkcjeKsztaltu_dS_fun[2][2]) + (funkcjeKsztaltu_dS_fun[3][2] * funkcjeKsztaltu_dS_fun[3][2]))*alfaFun *bokFun / 2;
	return matrixResult;
}
double** liczHpow3(double funkcjeKsztaltu_dS_fun[4][4], double alfaFun, double bokFun)
{
	double** matrixResult = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrixResult[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixResult[i][j] = 0;
		}
	}
	matrixResult[2][2] = ((funkcjeKsztaltu_dS_fun[4][2] * funkcjeKsztaltu_dS_fun[4][2]) + (funkcjeKsztaltu_dS_fun[5][2] * funkcjeKsztaltu_dS_fun[5][2]))*alfaFun *bokFun / 2;
	matrixResult[2][3] = ((funkcjeKsztaltu_dS_fun[4][2] * funkcjeKsztaltu_dS_fun[4][3]) + (funkcjeKsztaltu_dS_fun[5][2] * funkcjeKsztaltu_dS_fun[5][3]))*alfaFun *bokFun / 2;
	matrixResult[3][2] = ((funkcjeKsztaltu_dS_fun[4][3] * funkcjeKsztaltu_dS_fun[4][2]) + (funkcjeKsztaltu_dS_fun[5][3] * funkcjeKsztaltu_dS_fun[5][2]))*alfaFun *bokFun / 2;
	matrixResult[3][3] = ((funkcjeKsztaltu_dS_fun[4][3] * funkcjeKsztaltu_dS_fun[4][3]) + (funkcjeKsztaltu_dS_fun[5][3] * funkcjeKsztaltu_dS_fun[5][3]))*alfaFun *bokFun / 2;
	return matrixResult;
}
double** liczHpow4(double funkcjeKsztaltu_dS_fun[4][4], double alfaFun, double bokFun)
{
	double** matrixResult = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrixResult[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixResult[i][j] = 0;
		}
	}
	matrixResult[0][0] = ((funkcjeKsztaltu_dS_fun[6][0] * funkcjeKsztaltu_dS_fun[6][0]) + (funkcjeKsztaltu_dS_fun[7][0] * funkcjeKsztaltu_dS_fun[7][0]))*alfaFun *bokFun / 2;
	matrixResult[0][3] = ((funkcjeKsztaltu_dS_fun[6][0] * funkcjeKsztaltu_dS_fun[6][3]) + (funkcjeKsztaltu_dS_fun[7][0] * funkcjeKsztaltu_dS_fun[7][3]))*alfaFun *bokFun / 2;
	matrixResult[3][0] = ((funkcjeKsztaltu_dS_fun[6][3] * funkcjeKsztaltu_dS_fun[6][0]) + (funkcjeKsztaltu_dS_fun[7][3] * funkcjeKsztaltu_dS_fun[7][0]))*alfaFun *bokFun / 2;
	matrixResult[3][3] = ((funkcjeKsztaltu_dS_fun[6][3] * funkcjeKsztaltu_dS_fun[6][3]) + (funkcjeKsztaltu_dS_fun[7][3] * funkcjeKsztaltu_dS_fun[7][3]))*alfaFun *bokFun / 2;
	return matrixResult;
}
double* liczPpow1(double funkcjeKsztaltu_dS_fun[4][4], double ambientFun, double alfaFun, double bokFun)
{
	double* vectorResult = new double[4];
	for (int i = 0; i < 4; i++)
	{
		vectorResult[i] = 0;
	}
	vectorResult[0] = (-1)*(funkcjeKsztaltu_dS_fun[0][0] + funkcjeKsztaltu_dS_fun[1][0])*ambientFun*alfaFun *bokFun / 2;
	vectorResult[1] = (-1)*(funkcjeKsztaltu_dS_fun[0][1] + funkcjeKsztaltu_dS_fun[1][1])*ambientFun*alfaFun *bokFun / 2;
	return vectorResult;
}
double* liczPpow2(double funkcjeKsztaltu_dS_fun[4][4], double ambientFun, double alfaFun, double bokFun)
{
	double* vectorResult = new double[4];
	for (int i = 0; i < 4; i++)
	{
		vectorResult[i] = 0;
	}
	vectorResult[1] = (-1)*(funkcjeKsztaltu_dS_fun[2][1] + funkcjeKsztaltu_dS_fun[3][1])*ambientFun*alfaFun *bokFun / 2;
	vectorResult[2] = (-1)*(funkcjeKsztaltu_dS_fun[2][2] + funkcjeKsztaltu_dS_fun[3][2])*ambientFun*alfaFun *bokFun / 2;
	return vectorResult;
}
double* liczPpow3(double funkcjeKsztaltu_dS_fun[4][4], double ambientFun, double alfaFun, double bokFun)
{
	double* vectorResult = new double[4];
	for (int i = 0; i < 4; i++)
	{
		vectorResult[i] = 0;
	}
	vectorResult[2] = (-1)*(funkcjeKsztaltu_dS_fun[4][2] + funkcjeKsztaltu_dS_fun[5][2])*ambientFun*alfaFun *bokFun / 2;
	vectorResult[3] = (-1)*(funkcjeKsztaltu_dS_fun[4][3] + funkcjeKsztaltu_dS_fun[5][3])*ambientFun*alfaFun *bokFun / 2;
	return vectorResult;
}
double* liczPpow4(double funkcjeKsztaltu_dS_fun[4][4], double ambientFun, double alfaFun, double bokFun)
{
	double* vectorResult = new double[4];
	for (int i = 0; i < 4; i++)
	{
		vectorResult[i] = 0;
	}
	vectorResult[0] = (-1)*(funkcjeKsztaltu_dS_fun[6][0] + funkcjeKsztaltu_dS_fun[7][0])*ambientFun*alfaFun *bokFun / 2;
	vectorResult[3] = (-1)*(funkcjeKsztaltu_dS_fun[6][3] + funkcjeKsztaltu_dS_fun[7][3])*ambientFun*alfaFun *bokFun / 2;
	return vectorResult;
}
double* mnozMacierziWektor(double** macierz, double* wektor, int n)
{
	double* nowyWektor = new double[n];
	for (int i = 0; i < n; i++)
	{
		nowyWektor[i] = 0;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			nowyWektor[i] += macierz[i][j] * wektor[j];
		}
	}
	return nowyWektor;
}
//funkcja z eduinf roche przerobiona
bool gauss(int n, double ** AB, double * X)
{
	const double eps = 1e-12;
	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}
double liczMax(double** tab, int icfun, int wezly)
{
	double max = tab[icfun][0];
	
		for (int j = 0; j < wezly; j++)
		{
			if (tab[icfun][j] > max)
			{
				max = tab[icfun][j];
			}
		}
		return max;
}
double liczMin(double** tab, int icfun, int wezly)
{
	double min = tab[icfun][0];

	for (int j = 0; j < wezly; j++)
	{
		if (tab[icfun][j] < min)
		{
			min = tab[icfun][j];
		}
	}
	return min;
}