//
// Created by Xun Li on 2019-11-27.
//
#include <boost/algorithm/string.hpp>
#include <cmath>

#include "GenUtils.h"
#include "gda_data.h"

std::vector<std::vector<double> > gda_demean(const std::vector<std::vector<double> > &data) {
    std::vector<std::vector<double> > results(data.size());
    for (size_t i=0; i<data.size(); ++i) {
        results[i] = data[i];
        GenUtils::DeviationFromMean(results[i]);
    }
    return results;
}

std::vector<std::vector<double> > gda_standardize(const std::vector<std::vector<double> > &data) {
    std::vector<std::vector<double> > results(data.size());
    for (size_t i=0; i<data.size(); ++i) {
        results[i] = data[i];
        GenUtils::StandardizeData(results[i]);
    }
    return results;
}

std::vector<std::vector<double> > gda_standardize_mad(const std::vector<std::vector<double> > &data) {
    std::vector<std::vector<double> > results(data.size());
    for (size_t i=0; i<data.size(); ++i) {
        results[i] = data[i];
        GenUtils::MeanAbsoluteDeviation(results[i]);
    }
    return results;
}

std::vector<double> gda_naturalbreaks(int k, const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::NaturalBreaks(k, data, copy_undefs);
}

std::vector<double> gda_quantilebreaks(int k, const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::QuantileBreaks(k, data, copy_undefs);
}

std::vector<double> gda_hinge15breaks(const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::Hinge15Breaks(data, copy_undefs);
}

std::vector<double> gda_hinge30breaks(const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::Hinge30Breaks(data, copy_undefs);
}

std::vector<double> gda_percentilebreaks(const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::PercentileBreaks(data, copy_undefs);
}

std::vector<double> gda_stddevbreaks(const std::vector<double> &data, const std::vector<bool> &undefs) {
    std::vector<bool> copy_undefs = undefs;
    return GenUtils::StddevBreaks(data, copy_undefs);
}

void gda_transform_inplace(std::vector<double>& vals, const std::string& method)
{
    if (boost::iequals(method, "range_standardize")) {
        GenUtils::RangeStandardize(vals);
    } else if (boost::iequals(method, "range_adjust")) {
        GenUtils::RangeAdjust(vals);
    } else if (boost::iequals(method, "mad")) {
        GenUtils::MeanAbsoluteDeviation(vals);
    } else if (boost::iequals(method, "demean")) {
        GenUtils::DeviationFromMean(vals);
    } else {
        // z-standardization
        GenUtils::StandardizeData(vals);
    }
}

bool gda_rateStandardizeEB(const std::vector<double>& P,
                           const std::vector<double>& E,
                           std::vector<double>& results,
                           std::vector<bool>& undefined)
{
    int obs = (int)P.size();
    bool has_undef = false;
    double	sP=0.0, sE=0.0;
    double* p = new double[obs];
    int i = 0;

    // compute pi, the rate i, and the pop. rate b_hat
    for (i=0; i<obs; i++) {
        if (undefined[i]) {
            p[i] = 0;
            continue;
        }

        if (P[i] == 0.0) {
            undefined[i] = true;
            p[i] = 0;
        } else {
            sP += P[i];
            sE += E[i];
            p[i] = E[i] / P[i];
        }
    }

    if (sP == 0.0) {
        delete [] p;
        for (int i=0; i<obs; i++) {
            undefined[i] = true;
            results[i] = 0;
        }
        return has_undef;
    }

    const double b_hat = sE / sP;

    // compute a_hat, the variance
    double obs_valid = 0.0;
    double gamma=0.0;
    for (i=0; i< obs; i++) {
        if (!undefined[i]) {
            gamma += P[i] * ((p[i] - b_hat) * (p[i] - b_hat));
            obs_valid += 1;
        }
    }

    double a = (gamma / sP) - (b_hat / (sP / obs_valid));
    const double a_hat = a > 0 ? a : 0.0;

    for (i=0; i<obs; i++) {
        results[i] = 0.0;
        if (!undefined[i]) {
            const double se = P[i] > 0 ? sqrt(a_hat + b_hat/P[i]) : 0.0;
            results[i] = se > 0 ? (p[i] - b_hat) / se : 0.0;
        }
    }
    delete [] p;
    return !has_undef;
}
