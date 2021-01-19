//
// Created by Xun Li on 9/27/19.
//

#include <boost/algorithm/string.hpp>
#include "../GenUtils.h"
#include "DataUtils.h"
#include "cluster.h"
#include "redcap.h"
#include "redcap_wrapper.h"

using namespace SpanningTreeClustering;

redcap_wrapper::redcap_wrapper(unsigned int k,
        GeoDaWeight *w,
        const std::vector<std::vector<double> > &data,
        unsigned int redcap_method,
        const std::string &distance_method,
        const std::vector<double>& bound_vals,
        double min_bound,
        int rand_seed,
        int cpu_threads)
{
    if (w) {
        num_obs = w->num_obs;

        GalElement* gal = 0;
        gal = Gda::GetGalElement(w);
        if (gal) {
            // get control variable
            double *_bound_vals = 0;
            if ((int)bound_vals.size() == num_obs) {
                _bound_vals = new double[num_obs];
                for (int i = 0; i < num_obs; ++i) _bound_vals[i] = bound_vals[i];
            }
            // get distance matrix
            int n_cols = data.size();
            double** matrix = new double*[num_obs];
            int** mask = new int*[num_obs];
            for (int i=0; i<num_obs; ++i) {
                matrix[i] = new double[n_cols];
                mask[i] = new int[n_cols];
                for (int j=0; j<n_cols; ++j) mask[i][j] = 1.0;
            }
            for (int i=0; i<n_cols; ++i) {
                std::vector<double> vals = data[i];
                GenUtils::StandardizeData(vals);
                for (int r=0; r<num_obs; ++r) {
                    matrix[r][i] = vals[r];
                }
            }
            char dist = 'e';
            if (boost::iequals(distance_method, "manhattan")) dist = 'b';
            int transpose = 0; // row wise
            double* weight = new double[n_cols];
            for (int i=0; i<n_cols; ++i) weight[i] = 1.0;

            double** ragged_distances = distancematrix(num_obs, n_cols, matrix,  mask, weight, dist, transpose);
            double** distances = DataUtils::fullRaggedMatrix(ragged_distances, num_obs, num_obs);
            if (ragged_distances) {
                for (int i = 1; i < num_obs; i++) delete[] ragged_distances[i];
                delete[] ragged_distances;
            }

            // call redcap
            std::vector<bool> undefs(num_obs, false); // not used
            AbstractClusterFactory* redcap = 0;
            if (redcap_method == 0) {
                redcap = new FirstOrderSLKRedCap(num_obs, n_cols, distances, matrix, undefs,
                                                 gal, _bound_vals, min_bound, cpu_threads);
            } else if (redcap_method == 1) {
                redcap = new FullOrderCLKRedCap(num_obs, n_cols, distances, matrix, undefs,
                                                 gal, _bound_vals, min_bound, cpu_threads);
            } else if (redcap_method == 2) {
                redcap = new FullOrderALKRedCap(num_obs, n_cols, distances, matrix, undefs,
                                                gal, _bound_vals, min_bound, true, cpu_threads);
            } else if (redcap_method == 3) {
                redcap = new FullOrderSLKRedCap(num_obs, n_cols, distances, matrix, undefs,
                                                gal, _bound_vals, min_bound, cpu_threads);
            } else if (redcap_method == 4) {
                redcap = new FullOrderWardRedCap(num_obs, n_cols, distances, matrix, undefs,
                                                gal, _bound_vals, min_bound, cpu_threads);
            }
            if (redcap) {
                redcap->Partitioning(k);
                cluster_ids = redcap->GetRegions();
            }

            if (weight) delete[] weight;
            if (_bound_vals) delete[] _bound_vals;
            if (distances) {
                for (int i = 1; i < num_obs; i++) delete[] distances[i];
                delete[] distances;
            }
            if (matrix) {
                for (int i = 0; i < num_obs; ++i) delete[] matrix[i];
                delete[] matrix;
            }
        }
    }
}

redcap_wrapper::~redcap_wrapper() {

}

const std::vector<int> redcap_wrapper::GetFlatClusters() {
    return GenUtils::flat_2dclusters(num_obs, cluster_ids);
}

const std::vector<std::vector<int> > redcap_wrapper::GetClusters() {
    return cluster_ids;
}
