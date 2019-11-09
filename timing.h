/*
 * Copyright (c) 2019 by Yuhao Gu. All rights reserved.
 * E-Mail: yhgu2000@outlook.com
V2.0 */

#pragma once
#include <iostream>

using TestedFunc = void (*)();

struct TimingRecords
{
	TestedFunc _f;				// 被测函数的指针
	double* _data_p;			// 测试记录数组
	unsigned int _size;			// 测试记录条数
	const char* _name_ps;		// 被测函数名

public:
	TimingRecords();
	TimingRecords(TestedFunc f, unsigned int size, const char* name_ps = nullptr);
	TimingRecords(const TimingRecords& cpy);
	TimingRecords(TimingRecords&& mov) noexcept;
	~TimingRecords() noexcept;

public:
	TimingRecords& operator=(const TimingRecords& cpy);
	TimingRecords& operator=(TimingRecords&& mov) noexcept;
};

struct AnalysisReport
{
	TimingRecords* _records_p;	// 原计时记录的指针
	double _mean;				// 样本均值
	double _vari;				// 样本方差
	double _socm;				// 样本二阶中心矩
	struct {
		double _degree;			// 置信度
		double _lowerBound;		// 区间下界
		double _upperBound;		// 区间上界
	} _interval;				// 置信区间
};

struct ComparisionReport
{
	enum Results {
		kBalance = 0,
		kF1Faster = 1,
		kF2Faster = 2
	} _result;					// 比较结果
	AnalysisReport* _ar1_p;		// f1分析结果的指针
	AnalysisReport* _ar2_p;		// f2分析结果的指针
	double _level;				// 显著性水平
	double _faterBy;			// 快相对慢的速度比例
};

using CompareMethod = ComparisionReport::Results(*)(AnalysisReport&, AnalysisReport&, double);

// 单次计时
double timing(TestedFunc f);

// 多次计时
TimingRecords timing(
	TestedFunc f,
	unsigned int sampleNum,
	const char* name = nullptr,
	std::ostream* watcher_p = nullptr
);

// 分析计时结果
AnalysisReport analyze(TimingRecords& records, double degree = 0.95);

// 方法一：双正态总体检验
ComparisionReport::Results compare_method_1(AnalysisReport& f1, AnalysisReport& f2, double level);

// 方法二：配对检验
ComparisionReport::Results compare_method_2(AnalysisReport& f1, AnalysisReport& f2, double level);

// 比较两个分析结果
ComparisionReport compare(
	AnalysisReport& f1,
	AnalysisReport& f2,
	CompareMethod method_p,		// 比较方法
	double level = 0.05			// 显著性水平
);

// 打印计时分析报告
std::ostream& operator<<(std::ostream& out, const AnalysisReport& ar);

// 打印比较报告
std::ostream& operator<<(std::ostream& out, const ComparisionReport& cp);

// 自动分析函数并打印结果报告
void auto_timing(
	TestedFunc f,
	int sampleNum = 5,
	const char* fName = nullptr,
	std::ostream& out = std::cout
);

// 自动对比两个函数并打印结果报告
void auto_compare(
	TestedFunc f1,
	TestedFunc f2,
	int sampleNum = 5,
	const char* f1Name = nullptr,
	const char* f2Name = nullptr,
	std::ostream& out = std::cout
);
