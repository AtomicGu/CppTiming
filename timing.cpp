/*
 * Copyright (c) 2019 by Yuhao Gu. All rights reserved.
 * E-Mail: yhgu2000@outlook.com
V2.0 */

#include "timing.h"
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace std;

TimingRecords::TimingRecords() :
	_f(nullptr),
	_data_p(nullptr),
	_size(0),
	_name_ps(nullptr)
{}

TimingRecords::TimingRecords(TestedFunc f, unsigned int size, const char* name_ps) :
	_f(f),
	_data_p(new double[size]),
	_size(size),
	_name_ps(name_ps ? strcpy(new char[strlen(name_ps) + 1], name_ps) : nullptr)
{
}

TimingRecords::TimingRecords(const TimingRecords& cpy) :
	_f(cpy._f),
	_data_p(static_cast<double*>(memcpy(new double[cpy._size], cpy._data_p, sizeof(double)* cpy._size))),
	_size(cpy._size),
	_name_ps(cpy._name_ps ? strcpy(new char[strlen(cpy._name_ps) + 1], cpy._name_ps) : nullptr)
{
}

TimingRecords::TimingRecords(TimingRecords&& mov) noexcept :
	_f(mov._f),
	_data_p(mov._data_p),
	_size(mov._size),
	_name_ps(mov._name_ps)
{
	mov._data_p = nullptr;
	mov._name_ps = nullptr;
}

TimingRecords::~TimingRecords() noexcept
{
	if (_data_p)
		delete[] _data_p;
	if (_name_ps)
		delete[] _name_ps;
}

TimingRecords& TimingRecords::operator=(const TimingRecords& cpy)
{
	if (this == &cpy)
		return *this;
	this->~TimingRecords();
	_f = cpy._f;
	_data_p = static_cast<double*>(memcpy(new double[cpy._size], cpy._data_p, sizeof(double) * cpy._size));
	_size = cpy._size;
	_name_ps = cpy._name_ps ? strcpy(new char[strlen(cpy._name_ps) + 1], cpy._name_ps) : nullptr;
	return *this;
}

TimingRecords& TimingRecords::operator=(TimingRecords&& mov) noexcept
{
	this->~TimingRecords();
	_f = mov._f;
	_data_p = mov._data_p;
	mov._data_p = nullptr;
	_size = mov._size;
	_name_ps = mov._name_ps;
	mov._name_ps = nullptr;
	return *this;
}

// 正态分布分位点计算函数
double n_alpha(double a, double mean = 0, double vari = 1)
{
	if (a == 0.5) return 0;
	static const double b[] = {
		0.1570796288E1,
		0.3706987906E-1,
		-0.8364353589E-3,
		-0.2250947176E-3,
		0.6841218299E-5,
		0.5824238515E-5,
		-0.1045274970E-5,
		0.8360937017E-7,
		-0.3231081277E-8,
		0.3657763036E-10,
		0.6936233982E-12
	};
	double alpha = 0;
	bool opposite = false;
	if (a > 0.5)
	{
		opposite = true;
		a = 1 - a;
	}
	double y = -log(4 * a * (1 - a));
	double u = 0;
	for (int i = 0; i < 11; ++i)
	{
		u += b[i] * pow(y, i);
	}
	if (opposite)
		u = -sqrt(y * u);
	else
		u = sqrt(y * u);
	return u * sqrt(vari) + mean;
}

// t分布分位点计算函数（暂时不会算，只能固定alpha=0.05，先预留接口）
double t_alpha(double alpha, int n)
{
	static double table_005[] = {
		0,
		6.3138,2.9200,2.3534,2.1318,2.0150,
		1.9432,1.8946,1.8595,1.8331,1.8125,
		1.7959,1.7823,1.7709,1.7613,1.7531,
		1.7459,1.7396,1.7341,1.7291,1.7247,
		1.7207,1.7171,1.7139,1.7109,1.7081,
		1.7056,1.7033,1.7011,1.6991,1.6973,
		1.6955,1.6939,1.6924,1.6909,1.6896,
		1.6883,1.6871,1.6860,1.6849,1.6839,
		1.6829,1.6820,1.6811,1.6802,1.6794
	};
	if (n < 1 || n > 45)
		return 0;
	else
		return table_005[n];
}

ComparisionReport::Results compare_method_1(AnalysisReport& f1, AnalysisReport& f2, double level)
{
	AnalysisReport* faster, * slower;
	if (f1._mean < f2._mean)
	{
		faster = &f1;
		slower = &f2;
	}
	else
	{
		faster = &f2;
		slower = &f1;
	}
	double S = sqrt(
		((faster->_records_p->_size - 1) * faster->_vari + (slower->_records_p->_size - 1) * slower->_vari)
		/ (static_cast<double>(faster->_records_p->_size) + slower->_records_p->_size - 2)
	);
	double t = (faster->_mean - slower->_mean) / (S * sqrt(1 / faster->_records_p->_size + 1 / slower->_records_p->_size));
	level = 1 - level;
	if (t <= -t_alpha(level, faster->_records_p->_size + slower->_records_p->_size - 2))
	{
		if (faster == &f1)
			return ComparisionReport::Results::kF1Faster;
		else
			return ComparisionReport::Results::kF2Faster;
	}
	else
		return ComparisionReport::Results::kBalance;
}

ComparisionReport::Results compare_method_2(AnalysisReport& f1, AnalysisReport& f2, double level)
{
	const AnalysisReport* faster, * slower;
	if (f1._mean < f2._mean)
		faster = &f1, slower = &f2;
	else
		faster = &f2, slower = &f1;
	unsigned int sample_num = f1._records_p->_size;
	double* y_records = new double[sample_num];
	double mean = 0;
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		mean += (y_records[i] = faster->_records_p->_data_p[i] - slower->_records_p->_data_p[i]);
	}
	mean /= sample_num;
	double vari = 0;
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		vari += pow(y_records[i] - mean, 2);
	}
	double t = mean / (sqrt(vari) / sqrt(sample_num));
	level = 1 - level;
	if (t <= -t_alpha(level, f1._records_p->_size - 1))
	{
		if (faster == &f1)
			return ComparisionReport::Results::kF1Faster;
		else
			return ComparisionReport::Results::kF2Faster;
	}
	else
		return ComparisionReport::Results::kBalance;
}

double timing(TestedFunc f)
{
	using clk_period = chrono::high_resolution_clock::period;
	constexpr double kPeriod = static_cast<double>(clk_period::num) / clk_period::den;
	auto begin_time = chrono::high_resolution_clock::now();
	f();
	auto end_time = chrono::high_resolution_clock::now();
	return kPeriod * (end_time - begin_time).count();
}

TimingRecords timing(
	TestedFunc f,
	unsigned int sample_num,
	const char* name,
	std::ostream* watcher_p
)
{
	TimingRecords timing_records(f, sample_num, name);
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		if (watcher_p)*watcher_p << "Timing " << f << " : " << i << endl;
		timing_records._data_p[i] = timing(f);
	}
	return std::move(timing_records);
}

AnalysisReport analyze(TimingRecords& records, double degree)
{
	AnalysisReport report;
	report._records_p = &records;
	// 计算均值
	report._mean = 0;
	report._vari = 0;
	for (unsigned int i = 0; i < records._size; ++i)
		report._mean += records._data_p[i];
	report._mean /= records._size;
	// 计算方差和二阶中心距
	for (unsigned int i = 0; i < records._size; ++i)
		report._vari += pow(records._data_p[i] - report._mean, 2);
	report._socm = report._vari / records._size;
	report._vari /= records._size - 1.0;
	// 计算置信区间
	report._interval._degree = degree;
	report._interval._lowerBound = report._mean + sqrt(report._vari) * n_alpha(0.5 - degree / 2);
	report._interval._upperBound = report._mean + sqrt(report._vari) * n_alpha(0.5 + degree / 2);
	return std::move(report);
}

ComparisionReport compare(
	AnalysisReport& f1,
	AnalysisReport& f2,
	CompareMethod method_p,
	double level
)
{
	ComparisionReport cp;
	cp._ar1_p = &f1, cp._ar2_p = &f2;
	cp._level = level;
	cp._result = method_p(f1, f2, level);
	if (f1._mean < f2._mean)
		cp._faterBy = f2._mean / f1._mean - 1;
	else
		cp._faterBy = f1._mean / f2._mean - 1;
	return std::move(cp);
}

std::ostream& operator<<(std::ostream& out, const AnalysisReport& ar)
{
	out << "+----------- Timing Report -----------+\n"
		<< "|                                     |\n"
		<< "|   Function Name: " << left << setw(19);
	if (ar._records_p->_name_ps)
		out << ar._records_p->_name_ps;
	else
		out << ar._records_p->_f;
	out << "|\n"
		<< "|                                     |\n"
		<< "|   Sample Records:                   |\n";
	for (unsigned int i = 0; i < ar._records_p->_size; ++i)
	{
		out << setw(12) << "|" << left << setw(5) << i << left << setw(21) << ar._records_p->_data_p[i] << "|\n";
	}
	out << "|                                     |\n"
		<< "+-------------------------------------+\n"
		<< "|                                     |\n"
		<< "|   Analyzation:                      |\n"
		<< "|       mean:    " << setw(21) << ar._mean << "|\n"
		<< "|       vari:    " << setw(21) << ar._vari << "|\n"
		<< "|       socm:    " << setw(21) << ar._socm << "|\n"
		<< "|       ival(0.95):                   |\n"
		<< "|         [" << right << setw(11) << ar._interval._lowerBound << "," << left << setw(11) << ar._interval._upperBound << "]   |\n"
		<< "|                                     |\n"
		<< "+-------------------------------------+\n";
	return out;
}

std::ostream& operator<<(std::ostream& out, const ComparisionReport& cp)
{
	const auto& ar1 = *cp._ar1_p;
	const auto& ar2 = *cp._ar2_p;
	const auto& tr1 = *ar1._records_p;
	const auto& tr2 = *ar2._records_p;
	out << "+-------------- Comparision Report -------------+\n"
		<< "|                                               |\n"
		<< "|   f1 Name: " << left << setw(35);
	if (tr1._name_ps)
		out << tr1._name_ps;
	else
		out << tr1._f;
	out << "|\n"
		<< "|   f2 Name: " << setw(35);
	if (tr1._name_ps)
		out << tr1._name_ps;
	else
		out << tr1._f;
	out << "|\n"
		<< "|                                               |\n"
		<< "|   Sample Records:                             |\n"
		<< "|     index       f1              f2            |\n";
	for (unsigned int i = 0, j = 0; i < tr1._size || j < tr2._size; ++i, ++j)
	{
		out << setw(8) << "|" << setw(7) << i;
		out << setw(16);
		if (i < tr1._size)
			out << tr1._data_p[i];
		else
			out << "******";
		out << setw(16);
		if (i < tr2._size)
			out << tr2._data_p[i];
		else
			out << "******";
		out << " |\n";
	}
	out << "|                                               |\n"
		<< "+-----------------------------------------------+\n"
		<< "|                                               |\n"
		<< "|   Analysis:                                   |\n"
		<< "|     mean:    " << setw(16) << ar1._mean << setw(16) << ar2._mean << " |\n"
		<< "|     vari:    " << setw(16) << ar1._vari << setw(16) << ar2._vari << " |\n"
		<< "|     socm:    " << setw(16) << ar1._socm << setw(16) << ar2._socm << " |\n"
		<< "|     ival(0.95):                               |\n"
		<< "|              [" << setw(11) << ar1._interval._lowerBound << ",    [" << setw(11) << ar2._interval._lowerBound << ",   |\n"
		<< "|                " << right << setw(11) << ar1._interval._upperBound << "]     " << setw(11) << ar2._interval._upperBound << "]  |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+\n"
		<< "|                                               |\n"
		<< "|   Significance Level: " << left << setw(24) << cp._level << "|\n"
		<< "|   Faster By:" << right << setprecision(4) << setw(12) << cp._faterBy << "%                     |\n"
		<< "|   Time Saved:" << setprecision(4) << setw(12) << 1.0 / (cp._faterBy + 1) - 1 << "%                    |\n"
		<< "|   Conclusion:                                 |\n"
		<< "|          ******************************       |\n";
	switch (cp._result)
	{
	case ComparisionReport::kBalance:
		out << "|          * no significant differences *       |\n";
		break;
	case ComparisionReport::kF1Faster:
		out << "|          * f1 is significantly faster *       |\n";
		break;
	case ComparisionReport::kF2Faster:
		out << "|          * f2 is significantly faster *       |\n";
		break;
	}
	out << "|          ******************************       |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+\n";
	return out;
}

void auto_timing(
	TestedFunc f,
	int sampleNum,
	const char* fName,
	std::ostream& out
)
{
	TimingRecords tr = timing(f, sampleNum, fName, &out);
	AnalysisReport ar = analyze(tr);
	out << ar << endl;
}

void auto_compare(
	TestedFunc f1,
	TestedFunc f2,
	int sampleNum,
	const char* f1Name,
	const char* f2Name,
	ostream& out
)
{
	out << "Timing started\n"
		<< "Timing f1" << endl;
	TimingRecords tr1 = timing(f1, sampleNum, f1Name, &out);
	out << "Timing f2" << endl;
	TimingRecords tr2 = timing(f2, sampleNum, f2Name, &out);
	out << "Timing finished\n"
		<< "Analysis started" << endl;
	AnalysisReport fp1 = analyze(tr1);
	AnalysisReport fp2 = analyze(tr2);
	out << "Analysis finished\n"
		<< "Comparision started" << endl;
	ComparisionReport cr = compare(fp1, fp2, compare_method_1, 0.05);
	out << "Comparision finished\n"
		<< "*** All Procedure Finished ***\n" << std::endl;
	out << cr
		<< "|                                               |\n"
		<< "|   Method II :                                 |\n"
		<< "|          ******************************       |\n";
	switch (compare_method_2(fp1, fp2, 0.05))
	{
	case ComparisionReport::kBalance:
		out << "|          * no significant differences *       |\n";
		break;
	case ComparisionReport::kF1Faster:
		out << "|          * f1 is significantly faster *       |\n";
		break;
	case ComparisionReport::kF2Faster:
		out << "|          * f2 is significantly faster *       |\n";
		break;
	}
	out << "|          ******************************       |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+" << endl;
}
