#include "GwtWeight.h"
#include <algorithm>
#include <fstream>
#include <string.h>
#include <iomanip> 

using namespace std;

GwtElement::~GwtElement()
{
	if (data) delete [] data;
	nbrs = 0;
}

bool GwtElement::alloc(const int sz)
{
	if (data) delete [] data;
	if (sz > 0) {
		nbrs = 0;
		data = new GwtNeighbor[sz];
	}
	return !empty();
}

double GwtElement::SpatialLag(const std::vector<double>& x,
                              const bool std) const
{
	double lag= 0;
	int cnt = 0;
	for (cnt= Size() - 1; cnt >= 0; cnt--) {
		//lag += data[cnt].weight * x[ data[cnt].nbx ];
		lag += x[ data[cnt].nbx ];
	}
	if (std && Size() > 1) lag /= Size();
	return lag;
}

double GwtElement::SpatialLag(const double *x, const bool std) const  {
	double    lag= 0;
	int cnt = 0;
	for (cnt= Size() - 1; cnt >= 0; cnt--) {
		//lag += data[cnt].weight * x[ data[cnt].nbx ];
		lag += x[ data[cnt].nbx ];
	}
	if (std && Size() > 1) lag /= Size();
	return lag;
}

bool GwtElement::Check(long nbr_idx)
{
    for (long i=0; i<nbrs; ++i) {
        if (data[i].nbx == nbr_idx) return true;
    }
    return false;
}

std::vector<long> GwtElement::GetNbrs()
{
    std::vector<long> nbr_ids;
    for (long i=0; i<nbrs; ++i) {
        nbr_ids.push_back(data[i].nbx);
    }
    return nbr_ids;
}

std::vector<double> GwtElement::GetNbrWeights()
{
    std::vector<double> nbr_w;
    for (long i=0; i<nbrs; ++i) {
        nbr_w.push_back(data[i].weight);
    }
    return nbr_w;
}
////////////////////////////////////////////////////////////////////////////////
// GWTWeight
////////////////////////////////////////////////////////////////////////////////


void GwtWeight::Update(const std::vector<bool>& undefs)
{
    
}

const std::vector<long> GwtWeight::GetNeighbors(int obs_idx)
{
    return gwt[obs_idx].GetNbrs();
}

const std::vector<double> GwtWeight::GetNeighborWeights(int obs_idx)
{
    return gwt[obs_idx].GetNbrWeights();
}

bool GwtWeight::CheckNeighbor(int obs_idx, int nbr_idx)
{
    return gwt[obs_idx].Check(nbr_idx);
}

bool GwtWeight::HasIsolates(GwtElement *gwt, int num_obs)
{
	if (!gwt) return false;
	for (int i=0; i<num_obs; i++) {
        if (gwt[i].Size() <= 0)
            return true;
    }
	return false;
}

int GwtWeight::GetNbrSize(int obs_idx)
{
    return gwt[obs_idx].Size();
}

double GwtWeight::SpatialLag(int obs_idx, const std::vector<double> &data)
{
    return gwt[obs_idx].SpatialLag(data);
}

void GwtWeight::GetNbrStats()
{
    double empties = 0;
    for (int i=0; i<num_obs; i++) {
        if (gwt[i].Size() == 0)
        empties += 1;
    }
    sparsity = empties / (double)num_obs;
    // others
    int sum_nnbrs = 0;
    vector<int> nnbrs_array;
    std::map<int, int> e_dict;
    for (int i=0; i<num_obs; i++) {
        GwtNeighbor* nbrs = gwt[i].dt();
        int n_nbrs = 0;
        for (int j=0; j<gwt[i].Size();j++) {
            int nbr = nbrs[j].nbx;
            if (i != nbr) {
                n_nbrs++;
                e_dict[i] = nbr;
                e_dict[nbr] = i;
            }
        }
        sum_nnbrs += n_nbrs;
        if (i==0 || n_nbrs < min_nbrs) min_nbrs = n_nbrs;
        if (i==0 || n_nbrs > max_nbrs) max_nbrs = n_nbrs;
        nnbrs_array.push_back(n_nbrs);
    }
    //double n_edges = e_dict.size() / 2.0;
    density = 100.0* sum_nnbrs / (double)(num_obs * num_obs);
    
    if (num_obs > 0) mean_nbrs = sum_nnbrs / (double)num_obs;
    std::sort(nnbrs_array.begin(), nnbrs_array.end());
    if (num_obs % 2 ==0) {
        median_nbrs = (nnbrs_array[num_obs/2-1] + nnbrs_array[num_obs/2]) / 2.0;
    } else {
        median_nbrs = nnbrs_array[num_obs/2];
    }
}

bool
GwtWeight::Save(const char *ofname, const char *layer_name, const char *id_var_name, const std::vector<int> &id_vec) {
    std::ofstream out;
    out.open(ofname);
    if (!(out.is_open() && out.good())) return false;

    std::string out_layer_name = layer_name;
    const char *ptr = strstr(layer_name, " ");
    if (ptr != NULL) {
        // if layer_name contains an empty space, the layer name should be
        // braced with quotes "layer name"
        out_layer_name = "\"" + out_layer_name + "\"";
    }

    int num_obs = (int) id_vec.size();
    out << "0 " << num_obs << " " << layer_name;
    out << " " << id_var_name << endl;

    for (int i=0; i<num_obs; ++i) {
        for (long nbr=0; nbr<gwt[i].Size(); ++nbr) {
            const GwtNeighbor& current = gwt[i].elt(nbr);
            out << id_vec[i] << ' ' << id_vec[current.nbx];
            out << ' ' << setprecision(9) << setw(18)
                << current.weight << endl;
        }
    }
    return true;
}

bool GwtWeight::Save(const char *ofname, const char *layer_name, const char *id_var_name,
                     const std::vector<const char *> &id_vec) {
    std::ofstream out;
    out.open(ofname);
    if (!(out.is_open() && out.good())) return false;

    std::string out_layer_name = layer_name;
    const char *ptr = strstr(layer_name, " ");
    if (ptr != NULL) {
        // if layer_name contains an empty space, the layer name should be
        // braced with quotes "layer name"
        out_layer_name = "\"" + out_layer_name + "\"";
    }

    int num_obs = (int) id_vec.size();
    out << "0 " << num_obs << " " << layer_name;
    out << " " << id_var_name << endl;

    for (int i=0; i<num_obs; ++i) {
        for (long nbr=0; nbr<gwt[i].Size(); ++nbr) {
            const GwtNeighbor& current = gwt[i].elt(nbr);
            out << id_vec[i] << ' ' << id_vec[current.nbx];
            out << ' ' << setprecision(9) << setw(18)
                << current.weight << endl;
        }
    }
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// Gda Function
////////////////////////////////////////////////////////////////////////////////
GalElement* Gda::Gwt2Gal(const GwtElement* g, int num_obs)
{
    
    if (g == NULL) return 0;
    GalElement* gal = new GalElement[num_obs];
    for (int i=0; i<num_obs; ++i) {
        gal[i].SetSizeNbrs(g[i].Size());
        for (long nbr=0; nbr<g[i].Size(); ++nbr) {
            const GwtNeighbor& current = g[i].elt(nbr);
            gal[i].SetNbr(nbr, current.nbx);

        }
    }
    return gal;
}
