/*
 * Copyright (c) 2019 by Yuhao Gu. All rights reserved.
 * E-Mail: yhgu2000@outlook.com
V1.0.1 */

#include "timing.h"
#include <iomanip>
#include <cmath>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

using namespace std;

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
	{
		return 0;
	}
	else
	{
		return table_005[n];
	}
}

ComparisionReport::RESULTS compare_method_1(AnalysisReport& f1, AnalysisReport& f2, double level)
{
	AnalysisReport* faster, * slower;
	if (f1.mean < f2.mean)
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
		((faster->records_p->_size - 1) * faster->vari + (slower->records_p->_size - 1) * slower->vari)
		/ (static_cast<double>(faster->records_p->_size) + slower->records_p->_size - 2)
	);
	double t = (faster->mean - slower->mean) / (S * sqrt(1 / faster->records_p->_size + 1 / slower->records_p->_size));
	level = 1 - level;
	if (t <= -t_alpha(level, faster->records_p->_size + slower->records_p->_size - 2))
	{
		if (faster = &f1)
			return ComparisionReport::RESULTS::F1_FASTER;
		else
			return ComparisionReport::RESULTS::F2_FASTER;
	}
	else return ComparisionReport::RESULTS::BALANCE;
}

ComparisionReport::RESULTS compare_method_2(AnalysisReport& f1, AnalysisReport& f2, double level)
{
	const AnalysisReport* faster, * slower;
	if (f1.mean < f2.mean)
		faster = &f1, slower = &f2;
	else
		faster = &f2, slower = &f1;
	unsigned int sample_num = f1.records_p->_size;
	double* y_records = new double[sample_num];
	double mean = 0;
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		mean += (y_records[i] = faster->records_p->_p[i] - slower->records_p->_p[i]);
	}
	mean /= sample_num;
	double vari = 0;
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		vari += pow(y_records[i] - mean, 2);
	}
	double t = mean / (sqrt(vari) / sqrt(sample_num));
	level = 1 - level;
	if (t <= -t_alpha(level, f1.records_p->_size - 1))
	{
		if (faster == &f1)
			return ComparisionReport::RESULTS::F1_FASTER;
		else
			return ComparisionReport::RESULTS::F2_FASTER;
	}
	else return ComparisionReport::RESULTS::BALANCE;
}

double timing(pFunc f)
{
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	f();
	QueryPerformanceCounter(&t2);
	return (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart;
}

TimingRecords timing(pFunc f, unsigned int sample_num, std::ostream* watcher_p)
{
	TimingRecords timing_records(f, sample_num);
	for (unsigned int i = 0; i < sample_num; ++i)
	{
		if (watcher_p)* watcher_p << "Timing " << f << " : " << i << endl;
		timing_records._p[i] = timing(f);
	}
	return std::move(timing_records);
}

AnalysisReport analyze(TimingRecords& records, double degree)
{
	AnalysisReport report;
	report.records_p = &records;
	// 计算均值
	report.mean = 0;
	report.vari = 0;
	for (unsigned int i = 0; i < records._size; ++i)
		report.mean += records._p[i];
	report.mean /= records._size;
	// 计算方差和二阶中心距
	for (unsigned int i = 0; i < records._size; ++i)
		report.vari += pow(records._p[i] - report.mean, 2);
	report.socm = report.vari / records._size;
	report.vari /= records._size - 1.0;
	// 计算置信区间
	report.interval.degree = degree;
	report.interval.lower_bound = report.mean + sqrt(report.vari) * n_alpha(0.5 - degree / 2);
	report.interval.upper_bound = report.mean + sqrt(report.vari) * n_alpha(0.5 + degree / 2);
	return report;
}

ComparisionReport compare(AnalysisReport& f1, AnalysisReport& f2, pCompareMethod method_p, double level)
{
	ComparisionReport cp;
	cp.ar1_p = &f1, cp.ar2_p = &f2;
	cp.level = level;
	cp.result = method_p(f1, f2, level);
	if (f1.mean < f2.mean)
		cp.fater_by = f2.mean / f1.mean - 1;
	else
		cp.fater_by = f1.mean / f2.mean - 1;
	return std::move(cp);
}

std::ostream& operator<<(std::ostream& out, const AnalysisReport& tp)
{
	out << "+----------- Timing Report -----------+\n"
		<< "|                                     |\n"
		<< "|   Function: " << setw(8) << tp.records_p->_f << "              |\n"
		<< "|   Sample Records:                   |" << endl;
	for (unsigned int i = 0; i < tp.records_p->_size; ++i)
	{
		out << "|" << setw(12) << i << "    " << setw(12) << left << tp.records_p->_p[i] << right << setw(10) << "|" << endl;
	}
	out << "|                                     |\n"
		<< "+-------------------------------------+\n"
		<< "|                                     |\n"
		<< "|   Analyzation:                      |\n"
		<< "|       mean:" << setw(14) << tp.mean << setw(13) << "|\n"
		<< "|       vari:" << setw(14) << tp.vari << setw(13) << "|\n"
		<< "|       socm:" << setw(14) << tp.socm << setw(13) << "|\n"
		<< "|       ival(0.95):                   |\n"
		<< "|         [" << setw(11) << tp.interval.lower_bound << "," << setw(11) << tp.interval.upper_bound << "]   |\n"
		<< "|                                     |\n"
		<< "+-------------------------------------+" << endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const ComparisionReport& cp)
{
	out << "+-------------- Comparision Report -------------+\n"
		<< "|                                               |\n"
		<< "|   Sample Records:                             |\n"
		<< "|     index      f1(" << setw(8) << cp.ar1_p->records_p->_f << ")    f2(" << setw(8) << cp.ar2_p->records_p->_f << ")   |" << endl;
	for (unsigned int i = 0, j = 0; i < cp.ar1_p->records_p->_size || j < cp.ar2_p->records_p->_size; ++i, ++j)
	{
		out << "|" << setw(8) << i;
		out << left << "      " << setw(16);
		if (i < cp.ar1_p->records_p->_size)
			out << cp.ar1_p->records_p->_p[i];
		else
			out << "******";
		out << setw(16);
		if (i < cp.ar2_p->records_p->_size)
			out << cp.ar2_p->records_p->_p[i];
		else
			out << "******";
		out << right << setw(2) << "|" << endl;
	}
	out << "|                                               |\n"
		<< "+-----------------------------------------------+\n"
		<< "|                                               |\n"
		<< "|   Analysis:                                   |\n"
		<< "|     mean:" << "    " << setw(11) << left << cp.ar1_p->mean << "     " << setw(12) << cp.ar2_p->mean << right << setw(7) << "|\n"
		<< "|     vari:" << "    " << setw(11) << left << cp.ar1_p->vari << "     " << setw(12) << cp.ar2_p->vari << right << setw(7) << "|\n"
		<< "|     socm:" << "    " << setw(11) << left << cp.ar1_p->socm << "     " << setw(12) << cp.ar2_p->socm << right << setw(7) << "|\n"
		<< "|     ival(0.95):                               |\n"
		<< "|              [" << left << setw(11) << cp.ar1_p->interval.lower_bound << ",    [" << setw(11) << cp.ar2_p->interval.lower_bound << ",   |\n"
		<< "|                " << right << setw(11) << cp.ar1_p->interval.upper_bound << "]     " << setw(11) << cp.ar2_p->interval.upper_bound << "]  |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+\n"
		<< "|                                               |\n"
		<< "|   Significance Level : " << left << setw(10) << cp.level << "             |\n"
		<< "|   Faster By:" << right << setprecision(4) << setw(12) << cp.fater_by << "%                     |\n"
		<< "|   Time Saved:" << setprecision(4) << setw(12) << 1.0 / (cp.fater_by + 1) - 1 << "%                    |\n"
		<< "|   Conclusion:                                 |\n"
		<< "|          ******************************       |\n";
	switch (cp.result)
	{
	case ComparisionReport::BALANCE:
		out << "|          * no significant differences *       |\n";
		break;
	case ComparisionReport::F1_FASTER:
		out << "|          * f1 is significantly faster *       |\n";
		break;
	case ComparisionReport::F2_FASTER:
		out << "|          * f2 is significantly faster *       |\n";
		break;
	}
	out << "|          ******************************       |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+" << endl;
	return out;
}

void auto_timing(pFunc f, int sample_num, std::ostream& out)
{
	TimingRecords tr = timing(f, sample_num, &out);
	AnalysisReport ar = analyze(tr);
	out << ar << endl;
}

void auto_compare(pFunc f1, pFunc f2, int sample_num, ostream& out)
{
	out << "Timing started\n"
		<< "Timing f1" << endl;
	TimingRecords tr1 = timing(f1, sample_num, &out);
	out << "Timing f2" << endl;
	TimingRecords tr2 = timing(f2, sample_num, &out);
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
	case ComparisionReport::BALANCE:
		out << "|          * no significant differences *       |\n";
		break;
	case ComparisionReport::F1_FASTER:
		out << "|          * f1 is significantly faster *       |\n";
		break;
	case ComparisionReport::F2_FASTER:
		out << "|          * f2 is significantly faster *       |\n";
		break;
	}
	out << "|          ******************************       |\n"
		<< "|                                               |\n"
		<< "+-----------------------------------------------+" << endl;
}
