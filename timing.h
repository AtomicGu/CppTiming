/*
 * Copyright (c) 2019 by Yuhao Gu. All rights reserved.
 * E-Mail: yhgu2000@outlook.com
V1.0.0 */

#pragma once
#include <iostream>

typedef void(*pFunc)();

class TimingRecords
{
public:
	pFunc _f;					// 分析的函数指针
	double* _p;
	unsigned int _size;

public:
	TimingRecords() :_f(nullptr), _p(nullptr), _size(0) {};
	TimingRecords(pFunc f, unsigned int size) :
		_f(f),
		_p(new double[size]),
		_size(size)
	{}
	TimingRecords(const TimingRecords& cpy) :
		_f(cpy._f),
		_p(new double[cpy._size]),
		_size(cpy._size)
	{
		memcpy(_p, cpy._p, cpy._size);
	}
	TimingRecords(TimingRecords&& mov) noexcept :
		_f(mov._f),
		_p(mov._p),
		_size(mov._size)
	{
		mov._p = nullptr;
	}
	~TimingRecords() noexcept
	{
		if (_p)
			delete[] _p;
	}

public:
	TimingRecords& operator=(const TimingRecords& cpy)
	{
		this->~TimingRecords();
		_f = cpy._f;
		_p = new double[cpy._size];
		_size = cpy._size;
		memcpy(_p, cpy._p, cpy._size);
	}
	TimingRecords& operator=(TimingRecords&& mov) noexcept
	{
		this->~TimingRecords();
		_f = mov._f;
		_p = mov._p;
		_size = mov._size;
		mov._p = nullptr;
	}
};

struct AnalysisReport
{
	TimingRecords* records_p;	// 原计时记录的指针
	double mean;				// 样本均值
	double vari;				// 样本方差
	double socm;				// 样本二阶中心矩
	struct {
		double degree;			// 置信度
		double lower_bound;		// 区间下界
		double upper_bound;		// 区间上界
	} interval;					// 置信区间
};

struct ComparisionReport
{
	enum RESULTS {
		BALANCE = 0,
		F1_FASTER = 1,
		F2_FASTER = 2
	} result;				// 比较结果
	AnalysisReport* ar1_p;	// f1分析结果的指针
	AnalysisReport* ar2_p;	// f2分析结果的指针
	double level;			// 显著性水平
	double fater_by;		// 快相对慢的速度比例
};

typedef ComparisionReport::RESULTS(*pCompareMethod)(AnalysisReport&, AnalysisReport&, double);

// 单次计时
double timing(pFunc f);

// 多次计时
TimingRecords timing(pFunc f, unsigned int sample_num, std::ostream* watcher_p = nullptr);

// 分析计时结果
AnalysisReport analyze(TimingRecords& records, double degree = 0.95);

// 方法一：双正态总体检验
ComparisionReport::RESULTS compare_method_1(AnalysisReport& f1, AnalysisReport& f2, double level);

// 方法二：配对检验
ComparisionReport::RESULTS compare_method_2(AnalysisReport& f1, AnalysisReport& f2, double level);

// 比较两个分析结果
ComparisionReport compare(
	AnalysisReport& f1,
	AnalysisReport& f2,
	pCompareMethod method_p,	// 比较方法
	double level = 0.05			// 显著性水平
);

// 打印计时分析报告
std::ostream& operator<<(std::ostream& out, const AnalysisReport& tp);

// 打印比较报告
std::ostream& operator<<(std::ostream& out, const ComparisionReport& cp);

// 自动分析函数并打印结果
void auto_timing(
	pFunc f,
	int sample_num = 5,
	std::ostream & out = std::cout
);

// 自动对比两个函数并打印结果
void auto_compare(
	pFunc f1,
	pFunc f2,
	int sample_num = 5,
	std::ostream & out = std::cout
);
