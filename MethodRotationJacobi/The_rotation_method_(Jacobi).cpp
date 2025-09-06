#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <iomanip>

void Copy_Matrix(double** Original, double** Copy, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Copy[i][j] = Original[i][j];
		}
	}
}

void Matrix_Symmetry_Generation(double** Matrix, int n)
{
	std::mt19937 generator(static_cast<unsigned int>(std::time(nullptr)));
	std::uniform_real_distribution<double> distributation(1.0, 3.0);
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			Matrix[i][j] = std::round(distributation(generator) * 100) / 100.0;
			Matrix[j][i] = Matrix[i][j];
		}
		Matrix[i][i] = std::round(distributation(generator) * 100) / 100.0;
	}
}

bool Symmetry(double** Matrix, int n)
{
	double eps = 0.0000001;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (fabs(Matrix[i][j] - Matrix[j][i]) > eps)
			{
				return false;
			}
		}
	}
	return true;
}

void Print_Matrix(double** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << Matrix[i][j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void Scan_Matrix(double** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cin >> Matrix[i][j];
		}
	}
}

void Mult_Matrix(double** Result, double** Matrix1, double** Matrix2, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Result[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				Result[i][j] += (Matrix1[i][k] * Matrix2[k][j]);
			}
		}
	}
}

double T(double** Matrix_k, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				sum += (fabs(Matrix_k[i][j]) * fabs(Matrix_k[i][j]));
			}
		}
	}
	return sum;
}

void Method_Jacobi(double** Matrix_k, double** EigenVec, int n, double eps, int& iter)
{
	do
	{
		double** Matrix_k_1 = new double* [n];
		for (int i = 0; i < n; i++)
		{
			Matrix_k_1[i] = new double[n];
			for (int j = 0; j < n; j++)
			{
				Matrix_k_1[i][j] = Matrix_k[i][j];
			}
		}

		double M; int i_max = 0, j_max = 1; M = fabs(Matrix_k_1[0][1]);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if ((fabs(M) <= fabs(Matrix_k_1[i][j])) && (i != j))
				{
					M = fabs(Matrix_k_1[i][j]); i_max = i; j_max = j;
				}
			}
		}

		double Alpha = 0, P = 0;
		P = ((2 * Matrix_k_1[i_max][j_max]) / (Matrix_k_1[i_max][i_max] - Matrix_k_1[j_max][j_max]));
		Alpha = ((atan(P)) / 2);

		double** R = new double* [n];	double** R_1 = new double* [n]; double** Result = new double* [n];
		for (int i = 0; i < n; i++)
		{
			R[i] = new double[n];	R_1[i] = new double[n];	Result[i] = new double[n];
			for (int j = 0; j < n; j++)
			{
				if (i == j) { R[i][j] = 1; R_1[i][j] = 1; }
				else { R[i][j] = 0; R_1[i][j] = 0; }
				Result[i][j] = 0;
			}
		}

		R[i_max][i_max] = cos(Alpha); R[i_max][j_max] = -sin(Alpha);
		R[j_max][i_max] = sin(Alpha); R[j_max][j_max] = cos(Alpha);

		R_1[i_max][i_max] = cos(Alpha); R_1[i_max][j_max] = sin(Alpha);
		R_1[j_max][i_max] = -sin(Alpha); R_1[j_max][j_max] = cos(Alpha);

		//
		double** Result_EV = new double* [n];
		for (int i = 0; i < n; i++)
		{
			Result_EV[i] = new double[n];
			for (int j = 0; j < n; j++)
			{
				Result_EV[i][j] = 0;
			}
		}
		Mult_Matrix(Result_EV, EigenVec, R, n);
		Copy_Matrix(Result_EV, EigenVec, n);
		for (int i = 0; i < n; i++)
		{
			delete[]Result_EV[i];
		}
		delete[]Result_EV;
		//

		Mult_Matrix(Result, R_1, Matrix_k_1, n);
		Mult_Matrix(Matrix_k, Result, R, n);

		for (int i = 0; i < n; i++)
		{
			delete[]Matrix_k_1[i];	delete[]R[i];	delete[]R_1[i];	delete[]Result[i];
		}
		delete[]Matrix_k_1;	delete[]R;	delete[]R_1;	delete[]Result;

		iter++;
	} while (T(Matrix_k, n) > eps);
}

void Check(double** A, double** E, double* EigenValues, double** EigenVec, int n)
{
	int kol = 0;
	double** A_Copy = new double* [n];
	double* ResVec = new double[n];
	for (int i = 0; i < n; i++) A_Copy[i] = new double[n];
	Copy_Matrix(A, A_Copy, n);
	while (kol < n)
	{
		for (int i = 0; i < n; i++) ResVec[i] = 0;
		Copy_Matrix(A_Copy, A, n);
		for (int i = 0; i < n; i++)
		{
			A[i][i] -= (E[i][i] * EigenValues[kol]);
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				ResVec[i] += (A[i][j] * EigenVec[j][kol]);
			}
		}

		int LeftSum = 0;
		for (int i = 0; i < n; i++)
		{
			LeftSum += fabs(ResVec[i]);
		}
		std::cout << LeftSum << "=0\n";

		kol++;
	}
	Copy_Matrix(A_Copy, A, n);
}

void Test_1()
{
	int n;
	std::cout << "Введите размерность матрицы: ";
	std::cin >> n;

	double** Matrix = new double* [n];
	double** A = new double* [n];
	double** EigenVec = new double* [n];
	double** E = new double* [n];
	for (int i = 0; i < n; i++)
	{
		Matrix[i] = new double[n];
		EigenVec[i] = new double[n];
		E[i] = new double[n];
		A[i] = new double[n];
	}

	std::cout << "Заполните матрицу данными:\n";
	Scan_Matrix(Matrix, n);

	if (Symmetry(Matrix, n))
	{
		std::cout << "Исходная матрица, симметрична.\n";

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)	E[i][j] = 1;
				else E[i][j] = 0;
			}
		}

		Copy_Matrix(Matrix, A, n);
		Copy_Matrix(E, EigenVec, n);

		std::cout << "Задайте количество точностей: ";
		int count_eps;
		std::cin >> count_eps;
		double tenth = 0.1;
		double eps = 0.00001;
		int iterration;

		for (int h = 0; h < count_eps; h++)
		{
			eps *= tenth;
			iterration = 0;

			Method_Jacobi(Matrix, EigenVec, n, eps, iterration);

			std::cout << "\nСобственные значения:\n";
			double* EigenValues = new double[n];
			for (int i = 0; i < n; i++)
			{
				EigenValues[i] = Matrix[i][i];
				std::cout << "lyambda" << i + 1 << "=" << EigenValues[i] << "; ";
			}

			std::cout << "\n\nСобственные вектора (Матрица собственных векторов) :\n";
			Print_Matrix(EigenVec, n);

			std::cout << "Проверка:\n";
			Check(A, E, EigenValues, EigenVec, n);

			std::cout << "\nТочность: \n" << eps << "\n" << "Количество итераций: " << iterration << "\n";
			std::cout << "______________________________________________________________\n";
			Copy_Matrix(A, Matrix, n);
			Copy_Matrix(E, EigenVec, n);
		}
	}
	else
	{
		std::cout << "Метод не применим, матрица не симметрична.\n";
	}
}

void Test_2()
{
	int n;
	std::cout << "Введите размерность матрицы: ";
	std::cin >> n;

	double** Matrix_Gen = new double* [n];
	double** A = new double* [n];
	double** EigenVec = new double* [n];
	double** E = new double* [n];
	for (int i = 0; i < n; i++)
	{
		Matrix_Gen[i] = new double[n];
		EigenVec[i] = new double[n];
		E[i] = new double[n];
		A[i] = new double[n];
	}

	Matrix_Symmetry_Generation(Matrix_Gen, n);

	if (Symmetry(Matrix_Gen, n))
	{
		std::cout << "Исходная матрица, симметрична:\n";
		Print_Matrix(Matrix_Gen, n);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)	E[i][j] = 1;
				else E[i][j] = 0;
			}
		}

		Copy_Matrix(Matrix_Gen, A, n);
		Copy_Matrix(E, EigenVec, n);

		std::cout << "Задайте количество точностей: ";
		int count_eps;
		std::cin >> count_eps;
		double tenth = 0.1;
		double eps = 0.00001;
		int iterration;

		for (int h = 0; h < count_eps; h++)
		{
			eps *= tenth;
			iterration = 0;

			Method_Jacobi(Matrix_Gen, EigenVec, n, eps, iterration);

			std::cout << "\nСобственные значения:\n";
			double* EigenValues = new double[n];
			for (int i = 0; i < n; i++)
			{
				EigenValues[i] = Matrix_Gen[i][i];
				std::cout << "lmda" << i + 1 << "=" << EigenValues[i] << "; ";
			}

			std::cout << "\n\nСобственные вектора (Матрица собственных векторов) :\n";
			Print_Matrix(EigenVec, n);

			std::cout << "Проверка:\n";
			Check(A, E, EigenValues, EigenVec, n);

			std::cout << "\nТочность: \n" << eps << "\n" << "Количество итераций: " << iterration << "\n";
			std::cout << "______________________________________________________________\n";
			Copy_Matrix(A, Matrix_Gen, n);
			Copy_Matrix(E, EigenVec, n);
		}
	}
	else
	{
		std::cout << "Метод не применим, матрица не симметрична.";
	}
}

void main()
{
	setlocale(LC_ALL, "rus");

	int flag = 0;
	std::cout << "Выберите тест:\n1) Ввести самому (любую, либо из предложенных);\n2) Сгенерировать;\nНомер теста: ";
	std::cin >> flag;
	switch (flag)
	{
	case 0:
	{
		std::cout << "Error!\n";
		break;
	}
	case 1:
	{
		Test_1();
		break;
	}
	case 2:
	{
		Test_2();
		break;
	}
	}
}