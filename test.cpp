/*
 * Copyright (c) 2019 by Yuhao Gu. All rights reserved.
 * E-Mail: yhgu2000@outlook.com
V1.0 */

// 简单的准确率测验程序

#include <iostream>
#include <random>
#include <ctime>
#include "timing.h"

int foo1_time = 0;

void foo1()
{
	int b = 3;
	for (int i = 0; i < foo1_time; ++i)
	{
		b *= b;
	}
}

int foo2_time = 0;
void foo2()
{
	int b = 3;
	for (int i = 0; i < foo2_time; ++i)
	{
		b *= b;
	}
}

int test_time = 100;

int main()
{
	std::cout << "test time: ";
	std::cin >> test_time;
	srand(time(0));
	int method1Right = 0;
	int method2Right = 0;
	for (int i = 0; i < test_time; ++i)
	{
		foo1_time = rand() % 1000000;
		foo2_time = rand() % 1000000;
		auto tr1 = timing(foo1, 5);
		auto tr2 = timing(foo2, 5);
		auto tp1 = analyze(tr1);
		auto tp2 = analyze(tr2);
		int result_1 = compare_method_1(tp1, tp2, 0.95);
		int result_2 = compare_method_2(tp1, tp2, 0.95);
		if (foo1_time == foo2_time)
		{
			if (result_1 == ComparisionReport::Results::kBalance)
				++method1Right;
			if (result_2 == ComparisionReport::Results::kBalance)
				++method2Right;
		}
		else if (foo1_time > foo2_time)
		{
			if (result_1 == ComparisionReport::Results::kF2Faster)
				++method1Right;
			if (result_2 == ComparisionReport::Results::kF2Faster)
				++method2Right;
		}
		else
		{
			if (result_1 == ComparisionReport::Results::kF1Faster)
				++method1Right;
			if (result_2 == ComparisionReport::Results::kF1Faster)
				++method2Right;
		}
	}
	std::cout << "test " << test_time << " times" << std::endl;
	std::cout << "method1:\t" << method1Right << " right" << std::endl;
	std::cout << "method2:\t" << method2Right << " right" << std::endl;
}
