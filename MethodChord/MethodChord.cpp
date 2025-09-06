#include <iostream>
#include <stdio.h>
#include <math.h>
#include <utility>
#include <chrono>

double F(double x)
{
	return (sqrt(log10(x + 2)) - x);
}

double dF(double x)
{
	return ((1 / (2 * log(10) * (x + 2) * sqrt(log10(x + 2)))) - 1);
}

double ddF(double x)
{
	return (-(1 + 2 * log(x + 2)) / (4 * log(10) * log(10) * (x + 2) * (x + 2) * log10(x + 2) * sqrt(log10(x + 2))));
}

std::pair<double,double> MethodChord(double a, double b, double eps)
{
	auto start = std::chrono::high_resolution_clock::now();

	double t = 0.2;
	double m = INT_MAX, M = INT_MIN;
	double current = a;
	while (current <= (b + t))
	{
		if (M <= abs(dF(current))) M = abs(dF(current));
		if (m >= abs(dF(current))) m = abs(dF(current));
		current += t;
	}

	double Xn, Xn_1;
	int flag = 0;
	if ((ddF(b)*F(b)) > 0)
	{
		Xn_1 = a; flag = 1;
	}
	else if ((ddF(a) * F(a)) > 0)
	{
		Xn_1 = b; flag = 2;
	}
	else
	{
		std::cout << "error!" << std::endl;
	}

	int iteration = 0;

	do
	{
		iteration += 1;
		if (iteration != 1) { Xn_1 = Xn; };
		if (flag == 1)
		{
			Xn = Xn_1 - ((F(Xn_1)) / (F(b) - F(Xn_1))) * (b - Xn_1);
		}
		else if (flag == 2)
		{
			Xn = Xn_1 - ((F(Xn_1)) / (F(Xn_1) - F(a))) * (Xn_1 - a);
		}
		else
		{
			std::cout << "error!" << std::endl;
		}
	} while (eps <= (((M - m) / m) * (abs(Xn - Xn_1))));
	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	double seconds = duration.count();

	return { Xn,seconds };
}

void FindRoots(double eps)
{
	double start = -100;
	double end = 100;
	double step = 0.5;
	double ksi;
	double time;

	for (double a = start; a <= end; a += step)
	{
		double b = a + step;
		double Fa = F(a);
		double Fb = F(b);

		if ((Fa * Fb) < 0)
		{
			std::pair<double, double> p = MethodChord(a, b, eps);
			ksi = p.first;
			time = p.second;
			std::printf("Корень x=%.9f",ksi);
			std::cout << " найден на интервале [" << a << ";" << b << "]\n";
			std::printf("Время исполнения алгоритма для корня: %.9f секунд\n",time);
		}
	}
}

std::pair<double, double> dihit(double a, double b, double eps) {
	auto start = std::chrono::high_resolution_clock::now();
	double c;
	int k = 0;
	while (fabs(b - a) > eps) {
		k += 1;
		c = (a + b) / 2.0;
		if (F(a) * F(c) < 0) b = c;
		else if (F(b) * F(c) < 0) a = c;
		else {
			std::cout << "ERROR The root was not found";
			break;
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	double seconds = duration.count();
	return { c,seconds };
}

void main()
{
	setlocale(LC_ALL, "rus");

	std::cout << "Для проверки скорости сходимости в зависимости от заданной точности, введите 1, иначе 0.\n";
	int action;
	std::cin >> action;
	if (!action)
	{
		double eps;
		std::cout << "Введите точность:\n";
		std::cin >> eps;
		FindRoots(eps);
		std::cout << "Сравнение скорости сходимости от алгоритма представлено ниже, только введите параметры a,b\n";
		double a, b;
		std::cin >> a >> b;
		std::pair<double, double> p = dihit(a, b, eps);
		double ksi = p.first;
		double time = p.second;
		std::printf("Корень x=%.9f", ksi);
		std::cout << " найден на интервале [" << a << ";" << b << "]\n";
		std::printf("Время исполнения алгоритма для корня: %.9f секунд\n", time);
	}
	else
	{
		std::cout << "Введите количество точностей: ";
		int kol;
		std::cin >> kol;
		std::cout << "\n";
		std::cout << "Введите поочередно каждую точность:\n";
		double* E = new double[kol];
		for (int i = 0;i < kol;i++)
		{
			std::cout << "eps" << i+1 << "=";
			std::cin >> E[i];
		}

		std::cout << "Результаты представлены ниже\n";
		for (int i = 0; i < kol;i++)
		{
			std::cout << i + 1 << ") ";
			FindRoots(E[i]);
		}
	}

}