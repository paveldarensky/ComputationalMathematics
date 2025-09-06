#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <math.h>
#include <ctime>
#include <chrono>

void Print(double** Matrix_A, double* Vector_B, double* Vector_X, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%5.2f ", Matrix_A[i][j]);
		}
		printf("%5.2f \n", Vector_B[i]);
	}
	for (int i = 0; i < n; i++) printf("x%d=%.2f\n", i + 1, Vector_X[i]);
}

void Transformation(double** Matrix_A, double* Vector_B, int n)
{
	for (int i = 0; i < n; i++)
	{
		Vector_B[i] /= Matrix_A[i][i];
		for (int j = 0; j < n; j++)
		{
			if (i != j)	Matrix_A[i][j] /= (-Matrix_A[i][i]);
		}
		Matrix_A[i][i] = 0;
	}
}

bool CC(double** Matrix_A, int n)
{
	bool F1 = 1, F2 = 1, F3 = 1, F = 1;
	double one = 1;
	double MaxSum = 0;
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum += fabs(Matrix_A[i][j]);
		}
		MaxSum = std::max(MaxSum, sum);
	}
	F1 &= ((MaxSum < one) && (MaxSum > 0));
	F &= F1;

	MaxSum = 0;
	for (int j = 0; (j < n) && (F); j++)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += fabs(Matrix_A[i][j]);
		}
		MaxSum = std::max(MaxSum, sum);
	}
	F2 &= ((MaxSum < one) && (MaxSum > 0));
	F &= F2;

	double sum = 0;
	for (int j = 0; (j < n) && (F); j++)
	{
		for (int i = 0; i < n; i++)
		{
			sum += (fabs(Matrix_A[i][j]) * fabs(Matrix_A[i][j]));
		}
	}
	F3 &= ((sqrt(sum) < one) && (sum > 0));
	F &= F3;

	return F;
}

bool Condition(double* Vector_Xk, double* Vector_Xk_1, int n, double eps)
{
	bool F = 1;
	for (int i = 0; (i < n) && (F); i++)
	{
		F &= (((fabs(Vector_Xk[i] - Vector_Xk_1[i])) / (fabs(Vector_Xk[i]))) < eps);
	}
	return F;
}

void MSI(double** Matrix_A, double* Vector_B, double* Vector_Xk, int n, double eps)
{
	bool F;
	double* Vector_Xk_1 = new double[n];

	do
	{
		for (int i = 0; i < n; i++)	Vector_Xk_1[i] = Vector_Xk[i];

		for (int i = 0; i < n; i++)
		{
			Vector_Xk[i] = 0;
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				{
					Vector_Xk[i] += (Matrix_A[i][j] * Vector_Xk_1[j]);
				}
			}
			Vector_Xk[i] += Vector_B[i];
		}
		F = Condition(Vector_Xk, Vector_Xk_1, n, eps);
	} while (!F);
}

void RM(double** Matrix_A, double* Vector_C, double* Vector_Xk, int n, double eps)
{
	for (int i = 0; i < n; i++)
	{
		Vector_C[i] /= Matrix_A[i][i];
		for (int j = 0; j < n; j++)
		{
			if (i != j) Matrix_A[i][j] /= (-Matrix_A[i][i]);
		}
		Matrix_A[i][i] = -1;
	}

	double* R = new double[n];
	double sum_R = 0;
	for (int i = 0; i < n; i++) { R[i] = 1; sum_R += R[i]; }
	while (sum_R > eps)
	{
		for (int i = 0; i < n; i++)
		{
			double sum = 0;
			for (int j = 0; j < n; j++)
			{
				if (i != j) sum += (Matrix_A[i][j] * Vector_Xk[i]);
			}
			R[i] = Vector_C[i] - Vector_Xk[i] + sum;
		}

		double M = -1000; int ind = -1;
		for (int i = 0; i < n; i++) 
		{
			if (fabs(R[i]) >= M)
			{	
				M = R[i]; ind = i;
			}
		}

		for (int i = 0; i < n; i++) { if (i == ind) { Vector_Xk[i] += R[i]; } }
		sum_R = 0;
		for (int i = 0; i < n; i++) { sum_R += R[i]; }
	}
}

void Copy(double** Matrix_A, double* Vector_B, int n, double** Matrix_A_Copy, double* Vector_B_Copy)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Matrix_A_Copy[i][j] = Matrix_A[i][j];
		}
		Vector_B_Copy[i] = Vector_B[i];
	}
}

void Check(double** Matrix_A_Copy, double* Vector_B_Copy, double* Vector_X, int n)
{
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		double left_b = 0;
		for (int j = 0; j < n; j++)
		{
			left_b += (Matrix_A_Copy[i][j] * Vector_X[j]);
		}
		printf("%2.6f \"=\" %2.6f\n", left_b, Vector_B_Copy[i]);
	}
}

void Test_1()
{
	int n = 3;
	double** Matrix_A = new double* [n];
	double* Vector_B = new double[n];
	double* Vector_X = new double[n];
	double** Matrix_A_Copy = new double* [n];
	double* Vector_B_Copy = new double[n];
	for (int i = 0; i < n; i++) { Matrix_A[i] = new double[n]; Matrix_A_Copy[i] = new double[n]; }
	Matrix_A[0][0] = 4; Matrix_A[0][1] = 0.24; Matrix_A[0][2] = -0.08;
	Matrix_A[1][0] = 0.09; Matrix_A[1][1] = 3; Matrix_A[1][2] = -0.15;
	Matrix_A[2][0] = 0.04; Matrix_A[2][1] = -0.08; Matrix_A[2][2] = 4;
	Vector_B[0] = 8; Vector_B[1] = 9; Vector_B[2] = 20;
	Vector_X[0] = 0; Vector_X[1] = 0; Vector_X[2] = 0;
	
	Copy(Matrix_A, Vector_B, n, Matrix_A_Copy, Vector_B_Copy);
	Transformation(Matrix_A, Vector_B, n);
	Vector_X[0] = Vector_B[0]; Vector_X[1] = Vector_B[1]; Vector_X[2] = Vector_B[2];
	//Print(Matrix_A, Vector_B, Vector_X, n);
	if (CC(Matrix_A, n))
	{
		double eps;
		printf("Условие сходимости выполняется!\n");
		printf("Для исследования скорости зависимости от заданной точности введите 2, для одной точности (здесь же будет проверка) введите 1, для исследования скорости с другим алгоритмом введите 3: ");
		int flag = 0;
		scanf("%d", &flag);
		if (flag == 1)
		{
			printf("Задайте точность: ");
			std::cin >> eps;
			MSI(Matrix_A, Vector_B, Vector_X, n, eps);
			for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
			Check(Matrix_A_Copy, Vector_B_Copy, Vector_X, n);
		}
		else if (flag == 2)
		{
			printf("Введите количество точностей: ");
			int count = 0;
			std::cin >> count;
			eps = 0.1; double tenth = 0.1;
			for (int h = 0; h < count; h++)
			{
				for (int i = 0; i < n; i++) Vector_X[i] = Vector_B[i];
				auto start = std::chrono::high_resolution_clock::now();
				MSI(Matrix_A, Vector_B, Vector_X, n, eps);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> duration = end - start;
				double seconds = duration.count();
				std::cout << "_________________________________________\n";
				printf("Точность: %.9f\nВремя: %.9f\n", eps, seconds);
				for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
				std::cout << std::endl;
				eps *= tenth;
			}
		}
		else if (flag == 3)
		{
			printf("Задайте точность: ");
			std::cin >> eps;
			auto start1 = std::chrono::high_resolution_clock::now();
			MSI(Matrix_A, Vector_B, Vector_X, n, eps);
			auto end1 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> duration1 = end1 - start1;
			double seconds1 = duration1.count();
			std::cout << "_________________________________________\n";
			printf("Метод простой итерации\nВремя: %.9f\n", seconds1);
			for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
			for (int i = 0; i < n; i++) Vector_X[i] = 1;
			auto start2 = std::chrono::high_resolution_clock::now();
			RM(Matrix_A_Copy, Vector_B_Copy, Vector_X, n, eps);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> duration2 = end2 - start2;
			double seconds2 = duration2.count();
			std::cout << "_________________________________________\n";
			printf("Метод релаксации\nВремя: %.9f\n", seconds2);
			for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
		}
		else { std::cout << "error!"; }
	}
	else printf("Условие сходимости не выполняется!\n");
}

void Test_2()
{
	int n;
	printf("Введите количество неизвестных: ");
	std::cin >> n;
	printf("Заполните матрицу коэффициентов:\n");
	double** Matrix_A = new double* [n];
	double* Vector_B = new double[n];
	double* Vector_X = new double[n];
	double** Matrix_A_Copy = new double* [n];
	double* Vector_B_Copy = new double[n];
	for (int i = 0; i < n; i++)
	{
		Matrix_A[i] = new double[n]; Matrix_A_Copy[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			std::cin >> Matrix_A[i][j];
		}
	}
	printf("Заполните вектор свободных членов:\n");
	for (int i = 0; i < n; i++)	std::cin >> Vector_B[i];

	Copy(Matrix_A, Vector_B, n, Matrix_A_Copy, Vector_B_Copy);
	Transformation(Matrix_A, Vector_B, n);
	for (int i = 0; i < n; i++) Vector_X[i] = Vector_B[i];
	//Print(Matrix_A, Vector_B, Vector_X, n);
	if (CC(Matrix_A, n))
	{
		double eps;
		printf("Условие сходимости выполняется!\n");
		printf("Для исследования скорости зависимости от заданной точности введите 2, для одной точности (здесь же будет проверка) введите 1: ");
		int flag = 0;
		scanf("%d", &flag);
		if (flag == 1)
		{
			printf("Задайте точность: ");
			std::cin >> eps;
			MSI(Matrix_A, Vector_B, Vector_X, n, eps);
			for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
			Check(Matrix_A_Copy, Vector_B_Copy, Vector_X, n);
		}
		else if (flag == 2)
		{
			printf("Введите количество точностей: ");
			int count = 0;
			std::cin >> count;
			eps = 0.1; double tenth = 0.1;
			for (int h = 0; h < count; h++)
			{
				for (int i = 0; i < n; i++) Vector_X[i] = Vector_B[i];
				auto start = std::chrono::high_resolution_clock::now();
				MSI(Matrix_A, Vector_B, Vector_X, n, eps);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> duration = end - start;
				double seconds = duration.count();
				std::cout << "_________________________________________\n";
				printf("Точность: %.9f\nВремя: %.9f\n", eps, seconds);
				for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
				std::cout << std::endl;
				eps *= tenth;
			}
		}
		else { std::cout << "error!"; }
	}
	else printf("Условие сходимости не выполняется!\n");
}

void Test_3()
{
	int n;
	printf("Введите количество неизвестных: ");
	std::cin >> n;
	double** Matrix_A = new double* [n];
	double* Vector_B = new double[n];
	double* Vector_X = new double[n];
	double** Matrix_A_Copy = new double* [n];
	double* Vector_B_Copy = new double[n];
	std::srand(static_cast<unsigned int>(std::time(nullptr)));
	double start;
	double end;
	for (int i = 0; i < n; i++)
	{
		Matrix_A[i] = new double[n]; Matrix_A_Copy[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			if (j == i) { start = 60; end = 110; Matrix_A[i][j] = start + static_cast<double>(rand()) / RAND_MAX * (end - start); }
			else { start = -1; end = 1; Matrix_A[i][j] = start + static_cast<double>(rand()) / RAND_MAX * (end - start); }
		}
	}
	start = 10; end = 20;
	for (int i = 0; i < n; i++)	Vector_B[i] = start + static_cast<double>(rand()) / RAND_MAX * (end - start);
	Print(Matrix_A, Vector_B, Vector_X, n);
	Copy(Matrix_A, Vector_B, n, Matrix_A_Copy, Vector_B_Copy);
	Transformation(Matrix_A, Vector_B, n);
	for (int i = 0; i < n; i++) Vector_X[i] = Vector_B[i];
	//Print(Matrix_A, Vector_B, Vector_X, n);
	if (CC(Matrix_A, n))
	{
		double eps;
		printf("Условие сходимости выполняется!\n");
		printf("Для исследования скорости зависимости от заданной точности введите 2, для одной точности (здесь же будет проверка) введите 1: ");
		int flag = 0;
		scanf("%d", &flag);
		if (flag == 1)
		{
			printf("Задайте точность: ");
			std::cin >> eps;
			MSI(Matrix_A, Vector_B, Vector_X, n, eps);
			for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
			Check(Matrix_A_Copy, Vector_B_Copy, Vector_X, n);
		}
		else if (flag == 2)
		{
			printf("Введите количество точностей: ");
			int count = 0;
			std::cin >> count;
			eps = 0.1; double tenth = 0.1;
			for (int h = 0; h < count; h++)
			{
				for (int i = 0; i < n; i++) Vector_X[i] = Vector_B[i];
				auto start = std::chrono::high_resolution_clock::now();
				MSI(Matrix_A, Vector_B, Vector_X, n, eps);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> duration = end - start;
				double seconds = duration.count();
				std::cout << "_________________________________________\n";
				printf("Точность: %.9f\nВремя: %.9f\n", eps, seconds);
				for (int i = 0; i < n; i++) printf("x%d= %.6f\n", i + 1, Vector_X[i]);
				std::cout << std::endl;
				eps *= tenth;
			}
		}
		else { std::cout << "error!"; }
	}
	else printf("Условие сходимости не выполняется!\n");
}

void main()
{
	setlocale(LC_ALL, "rus");

	printf("Выберите тест:\n1 - СЛУ из Лабораторной работы №6;\n2 - СЛУ, заполненная Вами;\n3 - СЛУ, сгенерированная;\nВведите номер теста: ");
	int test;
	std::cin >> test;
	switch (test)
	{
	case 1:
	{
		printf("Тест %d\n", test);
		Test_1();
		break;
	}
	case 2:
	{
		printf("Тест %d\n", test);
		Test_2();
		break;
	}
	case 3:
	{
		printf("Тест %d\n", test);
		Test_3();
		break;
	}
	default:
		printf("Не выбран тест, либо неправильно введен номер теста!\n");
		break;
	}
}














