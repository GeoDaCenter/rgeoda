//
// Created by Xun Li on 2019-11-27.
//

#ifndef GEODA_GDA_DATA_H
#define GEODA_GDA_DATA_H

#include <vector>

// APIs for data processing

std::vector<std::vector<double> > gda_demean(const std::vector<std::vector<double> >& data);

std::vector<std::vector<double> > gda_standardize(const std::vector<std::vector<double> >& data);

std::vector<std::vector<double> > gda_standardize_mad(const std::vector<std::vector<double> >& data);

std::vector<double>  gda_naturalbreaks(int k, const std::vector<double>& data,
        const std::vector<bool>& undefs);

std::vector<double>  gda_quantilebreaks(int k, const std::vector<double>& data,
        const std::vector<bool>&  undefs);

std::vector<double>  gda_hinge15breaks(const std::vector<double>& data,
        const std::vector<bool>&  undefs);

std::vector<double>  gda_hinge30breaks(const std::vector<double>& data,
        const std::vector<bool>&  undefs);

std::vector<double>  gda_percentilebreaks(const std::vector<double>& data,
        const std::vector<bool>&  undefs);

std::vector<double>  gda_stddevbreaks( const std::vector<double>& data,
        const std::vector<bool>&  undefs);

void gda_transform_inplace(std::vector<double>& vals, const std::string& method);

bool gda_rateStandardizeEB(const std::vector<double>& P, const std::vector<double>& E,
                           std::vector<double>& results, std::vector<bool> &undefined);

#endif //GEODA_GDA_DATA_H
