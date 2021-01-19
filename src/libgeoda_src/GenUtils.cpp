/**
 * GeoDa TM, Copyright (C) 2011-2015 by Luc Anselin - all rights reserved
 *
 * This file is part of GeoDa.
 *
 * GeoDa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoDa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cfloat>
#include <set>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>

#include "rng.h"
#include "GdaConst.h"
#include "GenUtils.h"

using namespace std;

int StringUtils::utf8_strlen(const string& str)
{
    int c,i,ix,q;
    for (q=0, i=0, ix=str.length(); i < ix; i++, q++)
    {
        c = (unsigned char) str[i];
        if      (c>=0   && c<=127) i+=0;
        else if ((c & 0xE0) == 0xC0) i+=1;
        else if ((c & 0xF0) == 0xE0) i+=2;
        else if ((c & 0xF8) == 0xF0) i+=3;
        //else if (($c & 0xFC) == 0xF8) i+=4; // 111110bb //byte 5, unnecessary in 4 byte UTF-8
        //else if (($c & 0xFE) == 0xFC) i+=5; // 1111110b //byte 6, unnecessary in 4 byte UTF-8
        else return 0;//invalid utf8
    }
    return q;
}

void DbfFileUtils::SuggestDoubleParams(int length, int decimals,
                                       int* suggest_len, int* suggest_dec)
{
    // doubles have 52 bits for the mantissa, so we can allow at most
    // floor(log(2^52)) = 15 digits of precision.
    // We require that there length-2 >= decimals to allow for "x." . when
    // writing to disk, and when decimals = 15, require length >= 17 to
    // allow for "0." prefex. If length-2 == decimals, then negative numbers
    // are not allowed since there is not room for the "-0." prefix.
    //if (GdaConst::max_dbf_double_len < length) {
    if (35 < length) {
        length = 35;
    }
    if (length < 3) length = 3;
    if (decimals < 1) decimals = 1;
    if (decimals > 15) decimals = 15;
    if (length-2 < decimals) length = decimals + 2;

    *suggest_len = length;
    *suggest_dec = decimals;
}

double DbfFileUtils::GetMaxDouble(int length, int decimals,
                                  int* suggest_len, int* suggest_dec)
{
    // make sure that length and decimals have legal values
    SuggestDoubleParams(length, decimals, &length, &decimals);

    int len_inter = length - (1+decimals);
    //if (len_inter + decimals > 15) len_inter = 15-decimals;
    double r = 0;
    for (int i=0; i<len_inter+decimals; i++) r = r*10 + 9;
    for (int i=0; i<decimals; i++) r /= 10;

    if (suggest_len) *suggest_len = length;
    if (suggest_dec) *suggest_dec = decimals;
    return r;
}

std::string DbfFileUtils::GetMaxDoubleString(int length, int decimals)
{
    double x = GetMaxDouble(length, decimals, &length, &decimals);
    return GenUtils::doubleToString(x, decimals);
}

double DbfFileUtils::GetMinDouble(int length, int decimals,
                                  int* suggest_len, int* suggest_dec)
{
    SuggestDoubleParams(length, decimals, &length, &decimals);
    if (length-2 == decimals) return 0;
    if (suggest_len) *suggest_len = length;
    if (suggest_dec) *suggest_dec = decimals;
    return -DbfFileUtils::GetMaxDouble(length-1, decimals);
}

std::string DbfFileUtils::GetMinDoubleString(int length, int decimals)
{
    double x = GetMinDouble(length, decimals, &length, &decimals);
    if (length-2 == decimals) {
        std::string s("0.");
        for (int i=0; i<decimals; i++) s += "0";
        return s;
    }
    return GenUtils::doubleToString(x, decimals);
}

int DbfFileUtils::GetMaxInt(int length)
{
    // We want to allow the user to enter a string of
    // all 9s for the largest value reported.  So, we must
    // limit the length of the string to be floor(log(2^63)) = 18
    if (length < 1) return 0;
    if (length > 18) length = 18;
    int r=0;
    for (int i=0; i<length; i++) r = r*10 + 9;
    return r;
}

std::string DbfFileUtils::GetMaxIntString(int length)
{
    if (length < 19) {
        std::stringstream stream;
        stream << GetMaxInt(length);
        return stream.str();
    } else
        return "9223372036854775807"; // max value of int64
}

int DbfFileUtils::GetMinInt(int length)
{
    // This is generally the -GetMaxInt(length-1), because we must
    // allow one character for the minus sign unless the length
    // is greater than 18;
    if (length > 19) length = 19;
    return -GetMaxInt(length-1);
}

std::string DbfFileUtils::GetMinIntString(int length)
{
    if (length < 19) {
        std::stringstream stream;
        stream << GetMinInt(length);
        return stream.str();
    } else
        return "-9223372036854775808"; // min value of int64
}

std::string GenUtils::doubleToString(double val, int precision)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << val;
    return stream.str();
}

uint64_t Gda::ThomasWangHashUInt64(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

double Gda::ThomasWangHashDouble(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return 5.42101086242752217E-20 * key;
}

double Gda::ThomasWangDouble(uint64_t& key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return 5.42101086242752217E-20 * key;
}

double Gda::factorial(unsigned int n)
{
    double r=0;
    int i;
    for(i = n-1; i > 1; i--)
    r *= i;

    return r;
}

double Gda::nChoosek(unsigned int n, unsigned int k) {

    double r = 1;
    double s = 1;
    int i;
    int kk = k > n/2 ? k : n-k;

    for(i=n; i > kk; i--) r *= i;
    for(i=(n-kk); i>0; i--) s *= i;
    return r/s;
}

std::string Gda::CreateUUID(int nSize)
{
    if (nSize < 0 || nSize >= 38)
        nSize = 8;

    std::string letters = "abcdefghijklmnopqrstuvwxyz0123456789";

    Xoroshiro128Random rng;
    rng.SetSeed(4101842887655102017L);
    
    std::string uid;
    while ((int)uid.length() < nSize) {
        int iSecret = rng.nextLong() % letters.size();
        uid += letters[iSecret];
    }
    return uid;
}

/** Use with std::sort for sorting in ascending order */
bool Gda::dbl_int_pair_cmp_less(const dbl_int_pair_type& ind1,
								  const dbl_int_pair_type& ind2)
{
	return ind1.first < ind2.first;
}

/** Use with std::sort for sorting in descending order */
bool Gda::dbl_int_pair_cmp_greater(const dbl_int_pair_type& ind1,
									 const dbl_int_pair_type& ind2)
{
	return ind1.first > ind2.first;
}

/** Use with std::sort for sorting in ascending order */
bool Gda::dbl_int_pair_cmp_second_less(const dbl_int_pair_type& ind1,
										 const dbl_int_pair_type& ind2)
{
	return ind1.second < ind2.second;
}

/** Use with std::sort for sorting in descending order */
bool Gda::dbl_int_pair_cmp_second_greater(const dbl_int_pair_type& ind1,
											const dbl_int_pair_type& ind2)
{
	return ind1.second > ind2.second;
}


void
HingeStats::
CalculateHingeStats(const std::vector<Gda::dbl_int_pair_type>& data)
{
	num_obs = data.size();
	double N = num_obs;
	is_even_num_obs = (num_obs % 2) == 0;
	min_val = data[0].first;
	max_val = data[num_obs-1].first;
	Q2_ind = (N+1)/2.0 - 1;
	if (is_even_num_obs) {
		Q1_ind = (N+2)/4.0 - 1;
		Q3_ind = (3*N+2)/4.0 - 1;
	} else {
		Q1_ind = (N+3)/4.0 - 1;
		Q3_ind = (3*N+1)/4.0 - 1;
	}
	Q1 = (data[(int) floor(Q1_ind)].first +
		  data[(int) ceil(Q1_ind)].first)/2.0;
	Q2 = (data[(int) floor(Q2_ind)].first +
		  data[(int) ceil(Q2_ind)].first)/2.0;
	Q3 = (data[(int) floor(Q3_ind)].first +
		  data[(int) ceil(Q3_ind)].first)/2.0;
	IQR = Q3 - Q1;
	extreme_lower_val_15 = Q1 - 1.5*IQR;
	extreme_lower_val_30 = Q1 - 3.0*IQR;
	extreme_upper_val_15 = Q3 + 1.5*IQR;
	extreme_upper_val_30 = Q3 + 3.0*IQR;
	min_IQR_ind = -1;
	for (int i=0; i<num_obs; i++) {
		if (data[i].first < Q1) min_IQR_ind = i;
		else break;
	}
	if (min_IQR_ind < num_obs-1) min_IQR_ind++;
	max_IQR_ind = num_obs;
	for (int i=num_obs-1; i>=0; i--) {
		if (data[i].first > Q3) max_IQR_ind = i;
		else break;
	}
	if (max_IQR_ind > 0) max_IQR_ind--;
}

void
HingeStats::
CalculateHingeStats(const std::vector<Gda::dbl_int_pair_type>& data,
                    const std::vector<bool>& data_undef)
{
    num_obs = data.size();
    double N = 0.0;
    std::vector<double> data_valid;

    bool has_init = false;
    for (int i =0; i<num_obs; i++) {
        int obs_idx = data[i].second;
        if (!data_undef[obs_idx]) {
            double val = data[i].first;
            data_valid.push_back(val); // sorted
            if (!has_init) {
                min_val = val;
                max_val = val;
                has_init = true;
            }
            if (val < min_val)
                min_val = val;
            if (val > max_val)
                max_val = val;
        }
    }

    N = data_valid.size();

    is_even_num_obs = (data_valid.size() % 2) == 0;

    Q2_ind = (N+1)/2.0 - 1;
    if (is_even_num_obs) {
        Q1_ind = (N+2)/4.0 - 1;
        Q3_ind = (3*N+2)/4.0 - 1;
    } else {
        Q1_ind = (N+3)/4.0 - 1;
        Q3_ind = (3*N+1)/4.0 - 1;
    }

    if (N == 0 || N < Q3_ind) return;

    Q1 = (data_valid[(int) floor(Q1_ind)] + data_valid[(int) ceil(Q1_ind)])/2.0;
    Q2 = (data_valid[(int) floor(Q2_ind)] + data_valid[(int) ceil(Q2_ind)])/2.0;
    Q3 = (data_valid[(int) floor(Q3_ind)] + data_valid[(int) ceil(Q3_ind)])/2.0;

    IQR = Q3 - Q1;

    extreme_lower_val_15 = Q1 - 1.5*IQR;
    extreme_lower_val_30 = Q1 - 3.0*IQR;
    extreme_upper_val_15 = Q3 + 1.5*IQR;
    extreme_upper_val_30 = Q3 + 3.0*IQR;

    min_IQR_ind = -1;
    for (int i=0; i<num_obs; i++) {
        if (data[i].first < Q1) {
            min_IQR_ind = i;
        }
        else
            break;
    }
    if (min_IQR_ind < num_obs-1) {
        min_IQR_ind++;
    }
    max_IQR_ind = num_obs;

    for (int i=num_obs-1; i>=0; i--) {
        if (data[i].first > Q3) {
            max_IQR_ind = i;
        }
        else
            break;
    }
    if (max_IQR_ind > 0)
        max_IQR_ind--;
}

// Assume input v is sorted.  If not, can sort
// with std::sort(v.begin(), v.end())
// Testing: for v = {15, 20, 35, 40, 50},
// percentile(1, v) = 15, percentile(10, v) = 15, percentile(11) = 15.25
// percentile(50, v) = 35, percentile(89, v) = 49.5,
// percentile(90, v) = 50, percentile(99, v) = 50
double Gda::percentile(double x, const std::vector<double>& v)
{
	int N = v.size();
	double Nd = (double) N;
	double p_0 = (100.0/Nd) * (1.0-0.5);
	double p_Nm1 = (100.0/Nd) * (Nd-0.5);

    if (v.empty()) return 0;

	if (x <= p_0) return v[0];
	if (x >= p_Nm1) return v[N-1];

	for (int i=1; i<N; i++) {
		double p_i = (100.0/Nd) * ((((double) i)+1.0)-0.5);
		if (x == p_i) return v[i];
		if (x < p_i) {
			double p_im1 = (100.0/Nd) * ((((double) i))-0.5);
			return v[i-1] + Nd*((x-p_im1)/100.0)*(v[i]-v[i-1]);
		}
	}
	return v[N-1]; // execution should never get here
}

// Same assumptions as above
double Gda::percentile(double x, const Gda::dbl_int_pair_vec_type& v,
                       const std::vector<bool>& undefs)
{
    std::vector<double> valid_data;
    for (size_t i = 0; i<v.size(); i++ ) {
        double val = v[i].first;
        int ind = v[i].second;

        if (undefs[ind])
            continue;

        valid_data.push_back(val);
    }
    return percentile(x, valid_data);
}

// Same assumptions as above
double Gda::percentile(double x, const Gda::dbl_int_pair_vec_type& v)
{
	int N = v.size();
	double Nd = (double) N;
	double p_0 = (100.0/Nd) * (1.0-0.5);
	double p_Nm1 = (100.0/Nd) * (Nd-0.5);

	if (x <= p_0)
        return v[0].first;

	if (x >= p_Nm1)
        return v[N-1].first;

	for (int i=1; i<N; i++) {
		double p_i = (100.0/Nd) * ((((double) i)+1.0)-0.5);
		if (x == p_i)
            return v[i].first;
		if (x < p_i) {
			double p_im1 = (100.0/Nd) * ((((double) i))-0.5);
			return v[i-1].first + Nd*((x-p_im1)/100.0)*(v[i].first-v[i-1].first);
		}
	}
	return v[N-1].first; // execution should never get here
}


SampleStatistics::SampleStatistics()
	 : sample_size(0), min(0), max(0), mean(0),
    var_with_bessel(0), var_without_bessel(0),
    sd_with_bessel(0), sd_without_bessel(0)
{
}

SampleStatistics::SampleStatistics(const std::vector<double>& data)
	: sample_size(0), min(0), max(0), mean(0),
	var_with_bessel(0), var_without_bessel(0),
	sd_with_bessel(0), sd_without_bessel(0)
{
	CalculateFromSample(data);
}

SampleStatistics::SampleStatistics(const std::vector<double>& data,
                                   const std::vector<bool>& undefs)
	: sample_size(0), min(0), max(0), mean(0),
	var_with_bessel(0), var_without_bessel(0),
	sd_with_bessel(0), sd_without_bessel(0)
{
    std::vector<double> valid_data;
    for (size_t i=0; i<data.size(); i++) {
        if (undefs[i] == false)
            valid_data.push_back(data[i]);
    }
	CalculateFromSample(valid_data);
}

SampleStatistics::SampleStatistics(const std::vector<double>& data,
                                   const std::vector<bool>& undefs1,
                                   const std::vector<bool>& undefs2)
	: sample_size(0), min(0), max(0), mean(0),
	var_with_bessel(0), var_without_bessel(0),
	sd_with_bessel(0), sd_without_bessel(0)
{
    std::vector<double> valid_data;
    for (size_t i=0; i<data.size(); i++) {
        if (undefs1[i] || undefs2[i])
            continue;
        valid_data.push_back(data[i]);
    }
	CalculateFromSample(valid_data);
}

void SampleStatistics::CalculateFromSample(const std::vector<double>& data,
                                           const std::vector<bool>& undefs)
{
    std::vector<double> valid_data;
    for (size_t i=0; i<data.size(); i++) {
        if (undefs[i] == false)
            valid_data.push_back(data[i]);
    }
    CalculateFromSample(valid_data);
}
void SampleStatistics::CalculateFromSample(const std::vector<double>& data)
{
	sample_size = data.size();
	if (sample_size == 0) return;

	CalcMinMax(data, min, max);
	mean = CalcMean(data);

	double n = sample_size;
	double sum_squares = 0;
	for (int i=0, iend = data.size(); i<iend; i++) {
		sum_squares += data[i] * data[i];
	}

	var_without_bessel = (sum_squares/n) - (mean*mean);
	sd_without_bessel = sqrt(var_without_bessel);

	if (sample_size == 1) {
		var_with_bessel = var_without_bessel;
		sd_with_bessel = sd_without_bessel;
	} else {
		var_with_bessel = (n/(n-1)) * var_without_bessel;
		sd_with_bessel = sqrt(var_with_bessel);
	}
}

/** We assume that the data has been sorted in ascending order */
void
SampleStatistics::
CalculateFromSample(const std::vector<Gda::dbl_int_pair_type>& data_,
                    const std::vector<bool>& undefs)
{
    std::vector<double> data;
    for (int i=0, iend = data_.size(); i<iend; i++) {
        int id = data_[i].second;
        if (!undefs[id]) {
            data.push_back(data_[i].first);
        }
    }

	sample_size = data.size();
	if (sample_size == 0) return;

	min = data[0];
	max = data[sample_size-1];
	mean = CalcMean(data);

	double n = sample_size;
	double sum_squares = 0;
	for (int i=0, iend = data.size(); i<iend; i++) {
		sum_squares += data[i] * data[i];
	}

	var_without_bessel = (sum_squares/n) - (mean*mean);
	sd_without_bessel = sqrt(var_without_bessel);

	if (sample_size == 1) {
		var_with_bessel = var_without_bessel;
		sd_with_bessel = sd_without_bessel;
	} else {
		var_with_bessel = (n/(n-1)) * var_without_bessel;
		sd_with_bessel = sqrt(var_with_bessel);
	}
}

string SampleStatistics::ToString()
{
	ostringstream ss;
	ss << "sample_size = " << sample_size << endl;
	ss << "min = " << min << endl;
	ss << "max = " << max << endl;
	ss << "mean = " << mean << endl;
	ss << "var_with_bessel = " << var_with_bessel << endl;
	ss << "var_without_bessel = " << var_without_bessel << endl;
	ss << "sd_with_bessel = " << sd_with_bessel << endl;
	ss << "sd_without_bessel = " << sd_without_bessel << endl;
	return ss.str();
}

double SampleStatistics::CalcMin(const std::vector<double>& data)
{
	double min = std::numeric_limits<double>::max();
	for (int i=0, iend=data.size(); i<iend; i++) {
		if ( data[i] < min ) min = data[i];
	}
	return min;
}

double SampleStatistics::CalcMax(const std::vector<double>& data)
{
	double max = -std::numeric_limits<double>::max();
	for (int i=0, iend=data.size(); i<iend; i++) {
		if ( data[i] > max ) max = data[i];
	}
	return max;
}

void SampleStatistics::CalcMinMax(const std::vector<double>& data,
								  double& min, double& max)
{
	if (data.size() == 0) return;
	min = data[0];
	max = data[0];
	for (int i=1, iend=data.size(); i<iend; i++) {
		if ( data[i] < min ) {
			min = data[i];
		} else if ( data[i] > max ) {
			max = data[i];
		}
	}
}


double SampleStatistics::CalcMean(const std::vector<double>& data)
{
	if (data.size() == 0) return 0;
	double total = 0;
	for (int i=0, iend=data.size(); i<iend; i++) {
		total += data[i];
	}
	return total / (double) data.size();
}

double SampleStatistics::CalcMean(
							const std::vector<Gda::dbl_int_pair_type>& data)
{
	if (data.size() == 0) return 0;
	double total = 0;
	for (int i=0, iend=data.size(); i<iend; i++) {
		total += data[i].first;
	}
	return total / (double) data.size();
}

SimpleLinearRegression::SimpleLinearRegression(const std::vector<double>& X,
											   const std::vector<double>& Y,
											   double meanX, double meanY,
											   double varX, double varY)
	: n(0), covariance(0), correlation(0), alpha(0), beta(0), r_squared(0),
	std_err_of_estimate(0), std_err_of_beta(0), std_err_of_alpha(0),
	t_score_alpha(0), t_score_beta(0), p_value_alpha(0), p_value_beta(0),
	valid(false), valid_correlation(false), valid_std_err(false),
	error_sum_squares(0)
{
	CalculateRegression(X, Y, meanX, meanY, varX, varY);
}

SimpleLinearRegression::SimpleLinearRegression(const std::vector<double>& X,
											   const std::vector<double>& Y,
                                               const std::vector<bool>& X_undef,
                                               const std::vector<bool>& Y_undef,
											   double meanX, double meanY,
											   double varX, double varY)
	: n(0), covariance(0), correlation(0), alpha(0), beta(0), r_squared(0),
	std_err_of_estimate(0), std_err_of_beta(0), std_err_of_alpha(0),
	t_score_alpha(0), t_score_beta(0), p_value_alpha(0), p_value_beta(0),
	valid(false), valid_correlation(false), valid_std_err(false),
	error_sum_squares(0)
{

    std::vector<double> X_valid;
    std::vector<double> Y_valid;

    for (size_t i=0; i<X.size(); i++) {
        if (X_undef[i] || Y_undef[i])
            continue;

        X_valid.push_back(X[i]);
        Y_valid.push_back(Y[i]);
    }
	CalculateRegression(X_valid, Y_valid, meanX, meanY, varX, varY);
}

void SimpleLinearRegression::CalculateRegression(const std::vector<double>& X,
												 const std::vector<double>& Y,
												 double meanX, double meanY,
												 double varX, double varY)
{
    n = X.size();
	if (X.size() != Y.size() || X.size() < 2 )
        return;
	double expectXY = 0;
	for (int i=0, iend=X.size(); i<iend; i++) {
		expectXY += X[i]*Y[i];
	}
	expectXY /= (double) X.size();
	covariance = expectXY - meanX * meanY;
	if (varX > 4*DBL_MIN) {
		beta = covariance / varX;
		alpha = meanY - beta * meanX;
		valid = true;
	}
	double SS_tot = varY*Y.size();
	error_sum_squares = 0; // error_sum_squares = SS_err
	double err=0;
	for (int i=0, iend=Y.size(); i<iend; i++) {
		err = Y[i] - (alpha + beta * X[i]);
		error_sum_squares += err * err;
	}
	if (error_sum_squares < 16*DBL_MIN) {
		r_squared = 1;
	} else {
		r_squared = 1 - error_sum_squares / SS_tot;
	}

	if (Y.size()>2 && varX > 4*DBL_MIN) {
		// error_sum_squares/(n-k-1), k=1
		std_err_of_estimate = error_sum_squares/(Y.size()-2);
		std_err_of_estimate = sqrt(std_err_of_estimate);
		std_err_of_beta = std_err_of_estimate/sqrt(X.size()*varX);
		double sum_x_squared = 0;
		for (int i=0, iend=X.size(); i<iend; i++) {
			sum_x_squared += X[i] * X[i];
		}
		std_err_of_alpha = std_err_of_beta * sqrt(sum_x_squared / X.size());

		if (std_err_of_alpha >= 16*DBL_MIN) {
			t_score_alpha = alpha / std_err_of_alpha;
		} else {
			t_score_alpha = 100;
		}
		if (std_err_of_beta >= 16*DBL_MIN) {
			t_score_beta = beta	/ std_err_of_beta;
		} else {
			t_score_beta = 100;
		}
		p_value_alpha = TScoreTo2SidedPValue(t_score_alpha, X.size()-2);
		p_value_beta = TScoreTo2SidedPValue(t_score_beta, X.size()-2);

		valid_std_err = true;
	}

	double d = sqrt(varX)*sqrt(varY);
	if (d > 4*DBL_MIN) {
		correlation = covariance / d;
		valid_correlation = true;
	}
}

double SimpleLinearRegression::TScoreTo2SidedPValue(double tscore, int df)
{
	using namespace boost::math;
	students_t dist(df);
	// Cumulative Distribution Function evaluated at tscore
	if ( tscore >= 0) {
		return 2*(1.0-cdf(dist, tscore));
	} else {
		return 2*cdf(dist,tscore);
	}

}

string SimpleLinearRegression::ToString()
{
	ostringstream ss;
	ss << "covariance = " << covariance << endl;
	ss << "correlation = " << correlation << endl;
	ss << "alpha = " << alpha << endl;
	ss << "beta = " << beta << endl;
	ss << "r_squared = " << r_squared << endl;
	ss << "valid = " << (valid ? "true" : "false") << endl;
	ss << "valid_correlation = " << (valid_correlation ? "true" : "false")
		<< endl;
	ss << "error_sum_squares = " << error_sum_squares << endl;
	return ss.str();
}

AxisScale::AxisScale()
: data_min(0), data_max(0), scale_min(0), scale_max(0),
scale_range(0), tic_inc(0), p(0)
{
}

AxisScale::AxisScale(double data_min_s, double data_max_s, int ticks_s,
                     int lbl_precision_s, bool lbl_prec_fixed_point_s)
: data_min(0), data_max(0), scale_min(0), scale_max(0),
scale_range(0), tic_inc(0), lbl_precision(lbl_precision_s),
 lbl_prec_fixed_point(lbl_prec_fixed_point_s),  ticks(ticks_s), p(0)
{
	CalculateScale(data_min_s, data_max_s, ticks_s);
}

AxisScale::AxisScale(const AxisScale& s)
: data_min(s.data_min), data_max(s.data_max),
	scale_min(s.scale_min), scale_max(s.scale_max),
	scale_range(s.scale_range), tic_inc(s.tic_inc),
	lbl_precision(s.lbl_precision), lbl_prec_fixed_point(s.lbl_prec_fixed_point),
	ticks(s.ticks), p(s.p), tics(s.tics), tics_str(s.tics_str),tics_str_show(s.tics_str_show)
{
}

AxisScale& AxisScale::operator=(const AxisScale& s)
{
	data_min = s.data_min;
	data_max = s.data_max;
	scale_min = s.scale_min;
	scale_max = s.scale_max;
	scale_range = s.scale_range;
	tic_inc = s.tic_inc;
	p = s.p;
	tics = s.tics;
	tics_str = s.tics_str;
	tics_str_show = s.tics_str_show;
	ticks = s.ticks;
    lbl_precision = s.lbl_precision;
    lbl_prec_fixed_point = s.lbl_prec_fixed_point;
	return *this;
}

void AxisScale::CalculateScale(double data_min_s, double data_max_s,
							   const int ticks)
{
	if (data_min_s <= data_max_s) {
		data_min = data_min_s;
		data_max = data_max_s;
	} else {
		data_min = data_max_s;
		data_max = data_min_s;
	}

	double data_range = data_max - data_min;
	if ( data_range <= 2*DBL_MIN ) {
		scale_max = ceil((data_max + 0.05)*10)/10;
		scale_min = floor((data_min - 0.05)*10)/10;
		scale_range = scale_max - scale_min;
		p = 1;
		tic_inc = scale_range/2;
		tics.resize(3);
		tics_str.resize(3);
		tics[0] = scale_min;
		tics[1] = scale_min + tic_inc;
		tics[2] = scale_max;
	} else {
		p = floor(log10(data_range))-1;
		scale_max = ceil(data_max / pow((double)10,p)) * pow((double)10,p);
		scale_min = floor(data_min / pow((double)10,p)) * pow((double)10,p);
		scale_range = scale_max - scale_min;
		tic_inc = floor((scale_range / pow((double)10,p))/ticks)
			* pow((double)10,p);
		if (scale_min + tic_inc*(ticks+1) <= scale_max + 2*DBL_MIN) {
			tics.resize(ticks+2);
			tics_str.resize(ticks+2);
		} else {
			tics.resize(ticks+1);
			tics_str.resize(ticks+1);
		}
		for (int i=0, iend=tics.size(); i<iend; i++) {
			tics[i] = scale_min + i*tic_inc;
		}
	}
	tics_str_show.resize(tics_str.size());
	for (int i=0, iend=tics.size(); i<iend; i++) {
        tics_str[i] = GenUtils::DblToStr(tics[i], lbl_precision,
                                         lbl_prec_fixed_point);
		tics_str_show[i] = true;
	}
}

/** only display every other tic value */
void AxisScale::SkipEvenTics()
{
	for (size_t i=0; i<tics_str_show.size(); i++) tics_str_show[i] = (i%2 == 0);
}

void AxisScale::ShowAllTics()
{
	for (size_t i=0; i<tics_str_show.size(); i++) tics_str_show[i] = true;
}

string AxisScale::ToString()
{
	ostringstream ss;
	ss << "data_min = " << data_min << endl;
	ss << "data_max = " << data_max << endl;
	ss << "scale_min = " << scale_min << endl;
	ss << "scale_max = " << scale_max << endl;
	ss << "scale_range = " << scale_range << endl;
	ss << "p = " << p << endl;
	ss << "tic_inc = " << tic_inc << endl;
	for (int i=0, iend=tics.size(); i<iend; i++) {
		ss << "tics[" << i << "] = " << tics[i];
		ss << ",  tics_str[" << i << "] = " << tics_str[i] << endl;
	}
	ss << "Exiting AxisScale::CalculateScale" << endl;
	return ss.str();
}

std::string GenUtils::BoolToStr(bool b)
{
	return b ? "true" : "false";
}

bool GenUtils::StrToBool(const std::string& s)
{
	if (boost::iequals(s, "1")) return true;
    if (boost::iequals(s, "true")) return true;
	return false;
}

/** If input string has length < width, then prepends (or appends
 if pad_left=false) string with spaces so that total length is now width.
 If string length >= width, then returns original input string. */
std::string GenUtils::Pad(const std::string& s, int width, bool pad_left)
{
	if ((int)s.length() >= width) return s;
	int pad_len = width - s.length();
	std::stringstream output;
	if (!pad_left) output << s;
	for (int i=0; i<pad_len; i++) output << " ";
	if (pad_left) output << s;
	return output.str();
}

std::string GenUtils::PadTrim(const std::string& s, int width, bool pad_left)
{
    /*if (s.length() > width) {
        int trim_w = width - 2; //"xxx..xxx"
        int second_w = trim_w / 2;
        int first_w = trim_w - second_w;
        std::string tmp = s.SubString(0, first_w-2);
        tmp << ".." << s.SubString(s.length() - second_w -1, s.length()-1);
        return tmp;
    }
    int pad_len = width - s.length();
    std::string output;
    if (!pad_left) output << s;
    for (int i=0; i<pad_len; i++) output << " ";
    if (pad_left) output << s;
    return output;
     */
    return s;
}

std::string GenUtils::DblToStr(double x, int precision, bool fixed_point)
{
	std::stringstream ss;
    if (x < 10000000) {
        ss << std::fixed;
    }

     if (x == (int)x && fixed_point == false) {
         // The default should be that an integer is displayed as an integer
        ss << (int)x;
    } else {
        ss << std::setprecision(precision);
        ss << x;
    }

    //ss << std::setprecision(precision);
    //ss << x;
	return std::string(ss.str().c_str());
}

std::string GenUtils::IntToStr(int x, int precision)
{
    std::stringstream ss;

    if (x < 10000000) {
        ss << std::fixed;
    }
    ss << std::setprecision(precision);
    ss << x;

    return std::string(ss.str().c_str());
}


double GenUtils::Median(std::vector<double>& data)
{
    if (data.empty()) return 0;

    std::sort(data.begin(), data.end());

    int n = data.size();
    if (n % 2 == 1) return data[n/2];

    return 0.5 * (data[n/2 -1] + data[n/2]);
}

void GenUtils::DeviationFromMean(int nObs, double* data)
{
	if (nObs == 0) return;
	double sum = 0.0;
	for (int i=0, iend=nObs; i<iend; i++) sum += data[i];
	const double mean = sum / (double) nObs;
	for (int i=0, iend=nObs; i<iend; i++) data[i] -= mean;
}

void GenUtils::DeviationFromMean(int nObs, double* data, std::vector<bool>& undef)
{
	if (nObs == 0) return;

    int nValid = 0;
	double sum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i])
            continue;
        sum += data[i];
        nValid += 1;
    }
	const double mean = sum / (double) nValid;
    for (int i=0, iend=nObs; i<iend; i++) {
        data[i] -= mean;
    }
}

void GenUtils::DeviationFromMean(std::vector<double>& data)
{
	if (data.size() == 0) return;
	double sum = 0.0;
	for (int i=0, iend=data.size(); i<iend; i++) sum += data[i];
	const double mean = sum / (double) data.size();
	for (int i=0, iend=data.size(); i<iend; i++) data[i] -= mean;
}

void GenUtils::DeviationFromMean(std::vector<double>& data, std::vector<bool>& undef)
{
    if (data.size() == 0) return;
    double sum = 0.0;
    int n = 0;
    for (int i=0, iend=data.size(); i<iend; i++) {
        if (undef[i]) continue;
        sum += data[i];
        n++;
    }
    const double mean = sum / n;
    for (int i=0, iend=data.size(); i<iend; i++) {
        if (undef[i]) continue;
        data[i] -= mean;
    }
}

void GenUtils::MeanAbsoluteDeviation(int nObs, double* data)
{
    if (nObs == 0) return;
    double sum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) sum += data[i];
    const double mean = sum / (double) nObs;
    double mad = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        mad += std::abs(data[i] - mean);
    }
    mad = mad / nObs;
    if (mad == 0) return;
    for (int i=0, iend=nObs; i<iend; i++) {
        data[i] = (data[i] - mean) / mad;
    }
}

void GenUtils::MeanAbsoluteDeviation(int nObs, double* data,
                                     std::vector<bool>& undef)
{
    if (nObs == 0) return;
    double nValid = 0;
    double sum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i]) continue;
        sum += data[i];
        nValid += 1;
    }
    const double mean = sum / nValid;
    double mad = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i]) continue;
        mad += std::abs(data[i] - mean);
    }
    mad = mad / nValid;
    if (mad == 0) return;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i]) continue;
        data[i] = (data[i] - mean) / mad;
    }
}
void GenUtils::MeanAbsoluteDeviation(std::vector<double>& data)
{
	if (data.size() == 0) return;
	double sum = 0.0;
    double nn = data.size();
	for (int i=0, iend=data.size(); i<iend; i++) sum += data[i];
    const double mean = sum / nn;
    double mad = 0.0;
    for (int i=0, iend=data.size(); i<iend; i++) {
        mad += std::abs(data[i] - mean);
    }
    mad = mad / nn;
    if (mad == 0) return;
    for (int i=0, iend=data.size(); i<iend; i++) {
        data[i] = (data[i] - mean) / mad;
    }
}
void GenUtils::MeanAbsoluteDeviation(std::vector<double>& data,
                                     std::vector<bool>& undef)
{
    if (data.size() == 0) return;
    double sum = 0.0;
    double nValid = 0;
    for (int i=0, iend=data.size(); i<iend; i++) {
        if (undef[i]) continue;
        sum += data[i];
        nValid += 1;
    }
    const double mean = sum / nValid;
    double mad = 0.0;
    for (int i=0, iend=data.size(); i<iend; i++) {
        if (undef[i]) continue;
        mad += std::abs(data[i] - mean);
    }
    mad = mad / nValid;
    if (mad == 0) return;
    for (int i=0, iend=data.size(); i<iend; i++) {
        if (undef[i]) continue;
        data[i] = (data[i] - mean) / mad;
    }
}

void GenUtils::Transformation(int trans_type,
                              std::vector<std::vector<double> >& data,
                              std::vector<std::vector<bool> >& undefs)
{
    if (trans_type < 1) {
        return;
    }
    for (size_t i=0; i<data.size(); i++) {
        if (trans_type == 1) {
            // demean
            DeviationFromMean(data[i], undefs[i]);
        } else if (trans_type == 2) {
            // standarize (z)
            StandardizeData(data[i], undefs[i]);
        } else if (trans_type == 3) {
            // MAD
            MeanAbsoluteDeviation(data[i], undefs[i]);
        }
    }
}

double GenUtils::Correlation(std::vector<double>& x, std::vector<double>& y)
{
    int nObs = x.size();
    double sum_x = 0;
    double sum_y = 0;
    for (int i=0; i<nObs; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / nObs;
    double mean_y = sum_y / nObs;

    double ss_x = 0;
    double ss_y = 0;
    double ss_xy = 0;
    double d_x = 0, d_y = 0;
    for (int i=0; i<nObs; i++) {
        d_x = x[i] - mean_x;
        d_y = y[i] - mean_y;
        ss_x += d_x * d_x;
        ss_y += d_y * d_y;
        ss_xy += d_x * d_y;
    }

    double r = pow(ss_x * ss_y, 0.5);
    r = ss_xy / r;
    return r;
}

double GenUtils::Sum(std::vector<double>& data)
{
    double sum = 0;
    int nObs = data.size();
    for (int i=0; i<nObs; i++) sum += data[i];
    return sum;
}

double GenUtils::SumOfSquares(std::vector<double>& data)
{
    int nObs = data.size();
    if (nObs <= 1) return 0;
    GenUtils::DeviationFromMean(data);
    double ssum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) ssum += data[i] * data[i];
    return ssum;
}


double GenUtils::GetVariance(std::vector<double>& data)
{
    if (data.size() <= 1) return 0;
    GenUtils::DeviationFromMean(data);
    double ssum = 0.0;
    for (int i=0, iend=data.size(); i<iend; i++) ssum += data[i] * data[i];
    return ssum / data.size();
}

bool GenUtils::StandardizeData(int nObs, double* data)
{
	if (nObs <= 1) return false;
	GenUtils::DeviationFromMean(nObs, data);
	double ssum = 0.0;
	for (int i=0, iend=nObs; i<iend; i++) ssum += data[i] * data[i];
	const double sd = sqrt(ssum / (double) (nObs-1.0));
	if (sd == 0) return false;
	for (int i=0, iend=nObs; i<iend; i++) data[i] /= sd;
	return true;
}

bool GenUtils::StandardizeData(int nObs, double* data, std::vector<bool>& undef)
{
	if (nObs <= 1) return false;

    int nValid = 0;
    for (size_t i=0; i<undef.size(); i++) {
        if (!undef[i])
            nValid += 1;
    }

	GenUtils::DeviationFromMean(nObs, data, undef);
	double ssum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i])
            continue;
        ssum += data[i] * data[i];
    }
	const double sd = sqrt(ssum / (double) (nValid-1.0));
	if (sd == 0)
        return false;
    for (int i=0, iend=nObs; i<iend; i++) {
        data[i] /= sd;
    }
	return true;
}


bool GenUtils::StandardizeData(std::vector<double>& data)
{
	if (data.size() <= 1) return false;
	GenUtils::DeviationFromMean(data);
	double ssum = 0.0;
	for (int i=0, iend=data.size(); i<iend; i++) ssum += data[i] * data[i];
	const double sd = sqrt(ssum / (double) (data.size()-1.0));
	if (sd == 0) return false;
	for (int i=0, iend=data.size(); i<iend; i++) data[i] /= sd;
	return true;
}

bool GenUtils::StandardizeData(std::vector<double>& data, std::vector<bool>& undef)
{
    int nObs = data.size();
    if (nObs <= 1) return false;

    int nValid = 0;
    for (size_t i=0; i<undef.size(); i++) {
        if (!undef[i])
            nValid += 1;
    }

    GenUtils::DeviationFromMean(data, undef);
    double ssum = 0.0;
    for (int i=0, iend=nObs; i<iend; i++) {
        if (undef[i])
            continue;
        ssum += data[i] * data[i];
    }
    const double sd = sqrt(ssum / (double) (nValid-1.0));
    if (sd == 0)
        return false;
    for (int i=0, iend=nObs; i<iend; i++) {
        data[i] /= sd;
    }
    return true;
}


/*
 Reverse
 Changes the order of bytes in the presentation of a 4 byte number.
  */
int GenUtils::Reverse(const int &val)
{
	union {
		int v;
		char d[4];
	} chameleon;
	chameleon.v= val;
	char tmp = chameleon.d[0];
	chameleon.d[0] = chameleon.d[3];
	chameleon.d[3] = tmp;
	tmp = chameleon.d[1];
	chameleon.d[1] = chameleon.d[2];
	chameleon.d[2] = tmp;
	return chameleon.v;
}

long GenUtils::ReverseInt(const int &val)
{
	union {
		int v;
		char d[4];
	} chameleon;
	chameleon.v= val;
	char tmp = chameleon.d[0];
	chameleon.d[0] = chameleon.d[3];
	chameleon.d[3] = tmp;
	tmp = chameleon.d[1];
	chameleon.d[1] = chameleon.d[2];
	chameleon.d[2] = tmp;
	return chameleon.v;
}

void GenUtils::SkipTillNumber(std::istream &s)
{
	char ch;
	while (s >> ch) {
		if ((ch >= '0' && ch <= '9') || ch == '-' || ch == '+' || ch == '.')
			break;
	}
	if (s.good()) s.putback(ch);
}

// This is an implementation of ltoa
void GenUtils::longToString(const long d, char* Id, const int base)
{
	int i = 0;
	long j = d;
	//char rId[ GdaConst::ShpObjIdLen ];
	char rId[ 20 ];
	if (d == 0) {
		Id[0] = '0';
		Id[1] = '\0';
		return;
	}
	if (d < 0) j = -d;
	while (j != 0) {
		rId[i] = (j % base) + '0';
		j = j / base;
		i++;
	}
	j = i;
	if (d < 0) {
		Id[0] = '-';
		Id[i + 1] = '\0';
		while (i > 0) {
			Id[i] = rId[j - i];
			i--;
		}
		return;
	}

	Id[i] = '\0';
	while (i > 0) {
		Id[i - 1] = rId[j - i];
		i--;
	}
	return;
}



void GenUtils::strToInt64(const std::string& str, int *val)
{
	char buf[1024];
	strcpy( buf, (const char*)str.c_str());
	strToInt64(buf, val);
}

// Convert an ASCII string into a int (or long long)
void GenUtils::strToInt64(const char *str, int *val)
{
	int total = 0;
	bool minus = 0;

	while (isspace(*str)) str++;
	if (*str == '+') {
		str++;
	} else if (*str == '-') {
		minus = true;
		str++;
	}
	while (isdigit(*str)) {
		total *= 10;
		total += (*str++ - '0');
	}
	*val = minus ? -total : total;
}

bool GenUtils::validInt(const std::string& str)
{
	char buf[1024];
	strcpy( buf, (const char*)str.c_str());
	return validInt(buf);
}

// Checks that an ASCII string can be parsed to a valid integer.  At least
// one digit must been found.
bool GenUtils::validInt(const char* str)
{
	//LOG_MSG(std::string::Format("GenUtils::validInt(\"%s\"):", str));
	while (isspace(*str)) str++;
	if (*str == '+' || *str == '-') str++;
	const char* t = str;
	while (isdigit(*str)) str++;
	if (t == str) {
		// no digit found so return false
		//LOG_MSG("   no digit found");
		return false;
	}
	while (isspace(*str)) str++;
	// only return true if we are finally pointing at
	// the null terminating character.
	//LOG_MSG(std::string::Format("   final char is null: %d", *str == '\0'));
	return *str == '\0';
}

bool GenUtils::isEmptyOrSpaces(const std::string& str)
{
	char buf[1024];
	strcpy( buf, (const char*)str.c_str());
	return isEmptyOrSpaces(buf);
}

// returns true if the string is either empty
// or has only space characters
bool GenUtils::isEmptyOrSpaces(const char *str)
{
	while (isspace(*str)) str++;
	// if the first not-space char is not the end of the string,
	// return false.
	return *str == '\0';
}


std::string GenUtils::FindLongestSubString(const std::vector<std::string> strings,
										bool cs)
{
	using namespace std;
	int n=strings.size();
	if (n == 0) return "";
	vector<std::string> strs(strings);
	if (!cs) for (int i=0; i<n; i++)  boost::algorithm::to_lower(strs[i]);
	std::string ref_str = strs[0];
	for (int i=0; i<n; ++i) {
		if (strs[i].length() < ref_str.length()) ref_str = strs[i];
	}
	int len = ref_str.length();
	if (len == 0) return "";
	// iterate over all possible substrings in ref_str starting from first
	// position in ref_str, and starting with full ref_str.  Reduce length
	// of substring to search each iteration.
	for (int cur_len=len; cur_len > 0; --cur_len) {
		for (int cur_pos=0; cur_pos <= len-cur_len; ++cur_pos) {
			std::string ss = ref_str.substr(cur_pos, cur_len);
			bool all_match = true; // substring found everywhere currently
			for (int i=0; i<n && all_match; i++) {
				if (strs[i].find(ss) == std::string::npos) all_match = false;
			}
			if (all_match) {
				// common substring found.  Return unmodified (case-preserved)
				// substring from first string
				return strings[0].substr(strs[0].find(ss), cur_len);
			}
		}
	}
	return ""; // no substring match, return empty string.
}


bool GenUtils::less_vectors(const std::vector<int>& a,const std::vector<int>& b) {
    return a.size() > b.size();
}

const std::vector<int> GenUtils::flat_2dclusters(int n, std::vector<std::vector<int> > clusters) {
    vector<int> cluster_ids(n, 0);
    int ncluster = clusters.size();
    if (ncluster == 0)
        return cluster_ids;

    // sort result
    std::sort(clusters.begin(), clusters.end(), GenUtils::less_vectors);

    for (int i=0; i < ncluster; i++) {
        int c = i + 1;
        for (size_t j=0; j<clusters[i].size(); j++) {
            int idx = clusters[i][j];
            cluster_ids[idx] = c;
        }
    }
    return cluster_ids;
}

struct UniqueValElem {
    UniqueValElem(double v, int f, int l): val(v), first(f), last(l) {}
    double val; // value
    int first; // index of first occurrance
    int last; // index of last occurrance
};

/** clears uv_mapping and resizes as needed */
void create_unique_val_mapping(std::vector<UniqueValElem>& uv_mapping,
                               const std::vector<double>& v,
                               const std::vector<bool>& v_undef)
{
    uv_mapping.clear();
    int cur_ind = -1;

    for (int i=0, iend=v.size(); i<iend; i++) {
        if (v_undef[i])
            continue;
        if (uv_mapping.empty()) {
            cur_ind++;
            uv_mapping.push_back(UniqueValElem(v[i], i, i));
        } else {
            if (uv_mapping[cur_ind].val != v[i]) {
                uv_mapping[cur_ind].last = i-1;
                cur_ind++;
                uv_mapping.push_back(UniqueValElem(v[i], i, i));
            }
        }
    }
}

/** Assume that b.size() <= N-1 */
void pick_rand_breaks(std::vector<int>& b, int N, boost::uniform_01<boost::mt19937>& X)
{
    int num_breaks = b.size();
    if (num_breaks > N-1) return;

    std::set<int> s;
    while ((int)s.size() != num_breaks) s.insert(1 + (N-1)*X());
    int cnt=0;
    for (std::set<int>::iterator it=s.begin(); it != s.end(); it++) {
        b[cnt++] = *it;
    }
    std::sort(b.begin(), b.end());
}

/** Assume input b and v is sorted.  If not, can sort
 with std::sort(v.begin(), v.end())
 We assume that b and v are sorted in ascending order and are
 valid (ie, no break indicies out of range and all categories
 have at least one value.
 gssd is the global sum of squared differences from the mean */
double calc_gvf(const std::vector<int>& b, const std::vector<double>& v,
                double gssd)
{
    int N = v.size();
    int num_cats = b.size()+1;
    double tssd=0; // total sum of local sums of squared differences
    for (int i=0; i<num_cats; i++) {
        int s = (i == 0) ? 0 : b[i-1];
        int t = (i == num_cats-1) ? N : b[i];

        double m=0; // local mean
        double ssd=0; // local sum of squared differences (variance)
        for (int j=s; j<t; j++) m += v[j];
        m /= ((double) t-s);
        for (int j=s; j<t; j++) ssd += (v[j]-m)*(v[j]-m);
        tssd += ssd;
    }

    return 1-(tssd/gssd);
}

// translate unique value breaks into normal breaks given unique value mapping
void unique_to_normal_breaks(const std::vector<int>& u_val_breaks,
                             const std::vector<UniqueValElem>& u_val_mapping,
                             std::vector<int>& n_breaks)
{
    if (n_breaks.size() != u_val_breaks.size()) {
        n_breaks.resize(u_val_breaks.size());
    }
    for (int i=0, iend=u_val_breaks.size(); i<iend; i++) {
        n_breaks[i] = u_val_mapping[u_val_breaks[i]].first;
    }
}

std::vector<double>  GenUtils::NaturalBreaks(int num_cats, const vector<double>& data, vector<bool>& undef)
{

    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);

    std::vector<double> v(num_obs);
    std::vector<double> v_undef(num_obs);
    for (int i=0; i<num_obs; i++) {
        v[i] = var[i].first;
        int ind = var[i].second;
        v_undef[i] = undef[ind];
    }

    std::vector<UniqueValElem> uv_mapping;
    create_unique_val_mapping(uv_mapping, v, undef);

    int num_unique_vals = uv_mapping.size();
    int t_cats = std::min(num_unique_vals, num_cats);

    double mean = 0, max_val;
    int valid_obs = 0;
    for (int i=0; i<num_obs; i++) {
        double val = var[i].first;
        int ind = var[i].second;
        if (i==0 || val > max_val) max_val = val;
        if (undef[ind]) continue;
        mean += val;
        valid_obs += 1;
    }
    mean /= (double) valid_obs;

    double gssd = 0;
    for (int i=0; i<num_obs; i++) {
        double val = var[i].first;
        int ind = var[i].second;
        if (undef[ind]) continue;
        gssd += (val-mean)*(val-mean);
    }

    std::vector<int> rand_b(t_cats-1);
    std::vector<int> best_breaks(t_cats-1);
    std::vector<int> uv_rand_b(t_cats-1);

    double max_gvf_found = 0;
    int max_gvf_ind = 0;

    // for 5000 permutations, 2200 obs, and 4 time periods, slow enough
    // make sure permutations is such that this total is not exceeded.
    double c = 5000*2200*4;
    int perms = c / ((double) valid_obs);
    if (perms < 10) perms = 10;
    if (perms > 10000) perms = 10000;

    boost::mt19937 rng(GdaConst::gda_user_seed);
    boost::uniform_01<boost::mt19937> X(rng);

    for (int i=0; i<perms; i++) {
        pick_rand_breaks(uv_rand_b, num_unique_vals, X);
        // translate uv_rand_b into normal breaks
        unique_to_normal_breaks(uv_rand_b, uv_mapping, rand_b);
        double new_gvf = calc_gvf(rand_b, v, gssd);
        if (new_gvf > max_gvf_found) {
            max_gvf_found = new_gvf;
            max_gvf_ind = i;
            best_breaks = rand_b;
        }
    }

    std::vector<double> nat_breaks;
    nat_breaks.resize(best_breaks.size());
    for (int i=0, iend=best_breaks.size(); i<iend; i++) {
        nat_breaks[i] = var[best_breaks[i]].first;
    }

    return nat_breaks;
}

std::vector<double>  GenUtils::QuantileBreaks(int num_cats, const vector<double>& data, vector<bool>& undef) {

    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);


    std::vector<double> breaks(num_cats-1);
    for (int i=0, iend=breaks.size(); i<iend; i++) {
        breaks[i] = Gda::percentile(((i+1.0)*100.0)/((double) num_cats), var);
    }
    return breaks;
}

std::vector<double>  GenUtils::Hinge15Breaks(const vector<double>& data, vector<bool>& undef) {

    int num_cats = 6;
    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);


    std::vector<double> breaks(num_cats-1);

    HingeStats hinge_stats;
    hinge_stats.CalculateHingeStats(var);
    breaks[0] = hinge_stats.extreme_lower_val_15;
    breaks[1] = hinge_stats.Q1;
    breaks[2] = hinge_stats.Q2;
    breaks[3] = hinge_stats.Q3;
    breaks[4] = hinge_stats.extreme_upper_val_15;
    return breaks;
}

std::vector<double>  GenUtils::Hinge30Breaks(const vector<double>& data, vector<bool>& undef) {

    int num_cats = 6;
    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);


    std::vector<double> breaks(num_cats-1);

    HingeStats hinge_stats;
    hinge_stats.CalculateHingeStats(var);
    breaks[0] = hinge_stats.extreme_lower_val_30;
    breaks[1] = hinge_stats.Q1;
    breaks[2] = hinge_stats.Q2;
    breaks[3] = hinge_stats.Q3;
    breaks[4] = hinge_stats.extreme_upper_val_30;
    return breaks;
}

std::vector<double>  GenUtils::PercentileBreaks(const vector<double>& data, vector<bool>& undef) {

    int num_cats = 6;
    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);


    std::vector<double> breaks(num_cats-1);

    breaks[0] = Gda::percentile(1, var);
    breaks[1] = Gda::percentile(10, var);
    breaks[2] = Gda::percentile(50, var);
    breaks[3] = Gda::percentile(90, var);
    breaks[4] = Gda::percentile(99, var);
    return breaks;
}

std::vector<double>  GenUtils::StddevBreaks(const vector<double>& data, vector<bool>& undef) {

    int num_cats = 6;
    int num_obs = data.size();
    if (undef.size() == 0) {
        undef.resize(num_obs, false);
    }

    std::vector<std::pair<double, int> > var;
    for (int i=0; i<num_obs; ++i) {
        var.push_back(std::make_pair(data[i], i));
    }
    std::sort(var.begin(), var.end(), Gda::dbl_int_pair_cmp_less);


    std::vector<double> breaks(num_cats-1);

    std::vector<double> v(num_obs);
    SampleStatistics stats;
    for (int i=0; i<num_obs; i++) v[i] = var[i].first;
    stats.CalculateFromSample(v);
    breaks[0] = stats.mean - 2.0 * stats.sd_with_bessel;
    breaks[1] = stats.mean - 1.0 * stats.sd_with_bessel;
    breaks[2] = stats.mean;
    breaks[3] = stats.mean + 1.0 * stats.sd_with_bessel;
    breaks[4] = stats.mean + 2.0 * stats.sd_with_bessel;

    return breaks;
}
