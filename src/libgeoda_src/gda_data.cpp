//
// Created by Xun Li on 2019-11-27.
//

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

std::vector<double> gda_naturalbreaks(int k, const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::NaturalBreaks(k, data, copy_undefs);
}

std::vector<double> gda_quantilebreaks(int k, const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::QuantileBreaks(k, data, copy_undefs);
}

std::vector<double> gda_hinge15breaks(const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::Hinge15Breaks(data, copy_undefs);
}

std::vector<double> gda_hinge30breaks(const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::Hinge30Breaks(data, copy_undefs);
}

std::vector<double> gda_percentilebreaks(const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::PercentileBreaks(data, copy_undefs);
}

std::vector<double> gda_stddevbreaks(const std::vector<double> &data, const vector<bool> &undefs) {
    vector<bool> copy_undefs = undefs;
    return GenUtils::StddevBreaks(data, copy_undefs);
}
