#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <fstream>
#include <string.h>
#include <iomanip> 

#ifdef _WIN32
#if (_MSC_VER > 1900)
#include <functional>
#endif
#endif

#include "GalWeight.h"
#include "GwtWeight.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// GalElement
//
////////////////////////////////////////////////////////////////////////////////
GalElement::GalElement()
{
    is_nbrAvgW_empty = true;
}

bool GalElement::Check(long nbrIdx)
{
    if (nbrLookup.find(nbrIdx) != nbrLookup.end())
        return true;
    return false;
}

// return row standardized weights value
double GalElement::GetRW(int idx)
{
    if (is_nbrAvgW_empty) {
        size_t sz = nbr.size();
        nbrAvgW.resize(sz);
        double sumW = 0.0;
        
        for (size_t i=0; i<sz; i++)
            sumW += nbrWeight[i];
        
        for (size_t i=0; i<sz; i++) {
            nbrAvgW[i] = nbrWeight[i] / sumW;
        }
        is_nbrAvgW_empty = false;
    }
    
    if (nbrLookup.find(idx) != nbrLookup.end())
        return nbrAvgW[nbrLookup[idx]];
    return 0;
}

void GalElement::SetSizeNbrs(size_t	sz, bool is_gal)
{
	nbr.resize(sz);
    nbrWeight.resize(sz);
    if (!is_gal) {
        for (size_t i = 0; i < sz; i++) {
            nbrWeight[i] = 1.0;
        }
    }
}

// (which neighbor, what ID)
void GalElement::SetNbr(size_t pos, long n)
{
    if (pos < nbr.size()) {
        nbr[pos] = n;
        nbrLookup[n] = pos;
    }
    // this should be called by GAL created only
    if (pos < nbrWeight.size()) {
        nbrWeight[pos] = 1.0;
    }
}

// (which neighbor, what ID, what value)
void GalElement::SetNbr(size_t pos, long n, double w)
{
    if (pos < nbr.size()) {
        nbr[pos] = n;
        nbrLookup[n] = pos;
    } else {
        nbr.push_back(n);
        nbrLookup[n] = pos;
    }
    
    // this should be called by GWT-GAL 
    if (pos < nbrWeight.size()) {
        nbrWeight[pos] = w;
    } else {
        nbrWeight.push_back(w);
    }
}

// Update neighbor information on the fly using undefs information
// NOTE: this has to be used with a copy of weights (keep the original weights!)
void GalElement::Update(const std::vector<bool>& undefs)
{
    std::vector<int> undef_obj_positions;
   
    for (size_t i=0; i<nbr.size(); i++) {
        int obj_id = nbr[i];
        if (undefs[obj_id]) {
            int pos = nbrLookup[obj_id];
            undef_obj_positions.push_back(pos);
        }
    }
   
    if (undef_obj_positions.empty())
        return;
    
    // sort the positions in descending order, for removing from std::vector
	std::sort(undef_obj_positions.begin(),
              undef_obj_positions.end(), std::greater<int>());
   
    for (size_t i=0; i<undef_obj_positions.size(); i++) {
        size_t pos = undef_obj_positions[i];
        if (pos < nbr.size()) {
            nbrLookup.erase( nbr[pos] );
            nbr.erase( nbr.begin() + pos);
        }
        if (pos < nbrWeight.size()) {
            nbrWeight.erase( nbrWeight.begin() + pos);
        }
    }
}

void GalElement::SetNbrs(const GalElement& gal)
{
    size_t sz = gal.Size();
    nbr.resize(sz);
    nbrWeight.resize(sz);
    
    nbr = gal.GetNbrs();
    nbrLookup = gal.nbrLookup;
    nbrWeight = gal.GetNbrWeights();
    nbrLookup = gal.nbrLookup;
    nbrAvgW = gal.nbrAvgW;
}

const std::vector<long> & GalElement::GetNbrs() const
{
	return nbr;
}

const std::vector<double> & GalElement::GetNbrWeights() const
{
	return nbrWeight;
}

void GalElement::SortNbrs()
{
	std::sort(nbr.begin(), nbr.end(), std::greater<long>());
}

void GalElement::ReverseNbrs()
{
    std::reverse(nbr.begin(), nbr.end());
}

/** Compute spatial lag for a contiguity weights matrix.
 Automatically performs standardization of the result. */
double GalElement::SpatialLag(const std::vector<double>& x) const
{
	double lag = 0;
	size_t sz = Size();
   
    for (size_t i=0; i<sz; ++i) {
        lag += x[nbr[i]];
    }
    if (sz>1) lag /= (double) sz;
	
	return lag;
}

/** Compute spatial lag for a contiguity weights matrix.
 Automatically performs standardization of the result. */
double GalElement::SpatialLag(const double *x) const
{
	double lag = 0;
	size_t sz = Size();
    
    for (size_t i=0; i<sz; ++i) lag += x[nbr[i]];
    if (sz>1) lag /= (double) sz;

	return lag;
}

double GalElement::SpatialLag(const std::vector<double>& x,
							  const int* perm) const  
{
    // todo: this should also handle ReadGWtAsGAL like previous 2 functions
	double lag = 0;
	size_t sz = Size();
	for (size_t i=0; i<sz; ++i) lag += x[perm[nbr[i]]];
	if (sz>1) lag /= (double) sz;
	return lag;
}

////////////////////////////////////////////////////////////////////////////////
//
// GalWeight
//
////////////////////////////////////////////////////////////////////////////////
GalWeight::GalWeight(const GalWeight& gw)
: GeoDaWeight(gw)
{
	GalWeight::operator=(gw);
}

GalWeight& GalWeight::operator=(const GalWeight& gw)
{
	GeoDaWeight::operator=(gw);
	gal = new GalElement[num_obs];
    
    for (int i=0; i<num_obs; ++i) {
        gal[i].SetNbrs(gw.gal[i]);
    }
    
    this->num_obs = gw.num_obs;
    this->wflnm = gw.wflnm;
    this->id_field = gw.id_field;
    
	return *this;
}

void GalWeight::Update(const std::vector<bool>& undefs)
{
    for (int i=0; i<num_obs; ++i) {
        gal[i].Update(undefs);
    }

}

bool GalWeight::HasIsolates(GalElement *gal, int num_obs)
{
    if (!gal) {
        return false;
    }
	for (int i=0; i<num_obs; i++) {
        if (gal[i].Size() <= 0) {
            return true;
        }
    }
	return false;
}

int GalWeight::GetNbrSize(int obs_idx)
{
    return gal[obs_idx].Size();
}

double GalWeight::SpatialLag(int obs_idx, const std::vector<double> &data)
{
    return gal[obs_idx].SpatialLag(data);
}

void GalWeight::GetNbrStats()
{
    // sparsity
    double empties = 0;
    for (int i=0; i<num_obs; i++) {
        if (gal[i].Size() == 0)
        empties += 1;
    }
    sparsity = empties / (double)num_obs;
    
    // density
    // other
    int sum_nnbrs = 0;
    vector<int> nnbrs_array;
    std::map<int, int> e_dict;
    
    for (int i=0; i<num_obs; i++) {
        int n_nbrs = 0;
        const std::vector<long>& nbrs = gal[i].GetNbrs();
        for (size_t j=0; j<nbrs.size();j++) {
            int nbr = nbrs[j];
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
    density = 100.0 * sum_nnbrs / (double)(num_obs * num_obs);
    
    if (num_obs > 0) mean_nbrs = sum_nnbrs / (double)num_obs;
    std::sort(nnbrs_array.begin(), nnbrs_array.end());
    if (num_obs % 2 ==0) {
        median_nbrs = (nnbrs_array[num_obs/2-1] + nnbrs_array[num_obs/2]) / 2.0;
    } else {
        median_nbrs = nnbrs_array[num_obs/2];
    }
}

bool GalWeight::CheckNeighbor(int obs_idx, int nbr_idx)
{
    return gal[obs_idx].Check(nbr_idx);
}

const std::vector<long> GalWeight::GetNeighbors(int obs_idx)
{
    return gal[obs_idx].GetNbrs();
}

const std::vector<double> GalWeight::GetNeighborWeights(int obs_idx)
{
    return gal[obs_idx].GetNbrWeights();
}

bool GalWeight::Save(const char* ofname,
                     const char* layer_name,
                     const char* id_var_name,
                     const std::vector<int>& id_vec)
{
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
    size_t num_obs = (int) id_vec.size();
    out << "0 " << num_obs << " " << layer_name;
    out << " " << id_var_name << endl;

    for (size_t i=0; i<num_obs; ++i) {
        out << id_vec[i];
        out << " " << gal[i].Size() << endl;
        for (int cp=gal[i].Size(); --cp >= 0;) {
            out << id_vec[gal[i][cp]];
            if (cp > 0)
                out << " ";
        }
        out << endl;
    }
    return true;
}

bool GalWeight::Save(const char* ofname,
                     const char* layer_name,
                     const char* id_var_name,
                     const std::vector<const char*>& id_vec)
{
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

    size_t num_obs = (int) id_vec.size();
    out << "0 " << num_obs << " " << layer_name;
    out << " " << id_var_name << endl;

    for (size_t i=0; i<num_obs; ++i) {
        out << id_vec[i];
        out << " " << gal[i].Size() << endl;
        for (int cp=gal[i].Size(); --cp >= 0;) {
            out << id_vec[gal[i][cp]];
            if (cp > 0)
                out << " ";
        }
        out << endl;
    }
    return true;
}


/** Add higher order neighbors up to (and including) distance.
 If cummulative true, then include lower orders as well.  Otherwise,
 only include elements on frontier. */
void Gda::MakeHigherOrdContiguity(size_t distance, size_t obs,
                                  GalElement* W,
                                  bool cummulative)
{	
	if (obs < 1 || distance <=1) return;
	vector<vector<long> > X(obs);
	for (size_t i=0; i<obs; ++i) {
		vector<set<long> > n_at_d(distance+1);
		n_at_d[0].insert(i);
		for (size_t j=0, sz=W[i].Size(); j<sz; ++j) {
			n_at_d[1].insert(W[i][j]);
		}
		for (size_t d=2; d<=distance; ++d) {
			for (set<long>::const_iterator it=n_at_d[d-1].begin();
					 it!=n_at_d[d-1].end(); ++it)
			{
				for (size_t j=0, sz=W[*it].Size(); j<sz; ++j) {
					long nbr = W[*it][j];
					if (n_at_d[d-1].find(nbr) == n_at_d[d-1].end() &&
							n_at_d[d-2].find(nbr) == n_at_d[d-2].end()) {
						n_at_d[d].insert(nbr);
					}
				}
			}
		}
		size_t sz_Xi = 0;
		for (size_t d=(cummulative ? 1 : distance); d<=distance; ++d) {
			sz_Xi += n_at_d[d].size();
		}
		X[i].resize(sz_Xi);
		size_t cnt=0;
		for (size_t d=(cummulative ? 1 : distance); d<=distance; ++d) {
			for (set<long>::const_iterator it=n_at_d[d].begin();
					 it!=n_at_d[d].end(); ++it) { X[i][cnt++] = *it; }
		}
		sort(X[i].begin(), X[i].end(), greater<long>());
	}
	for (size_t i=0; i<obs; ++i) {
		W[i].SetSizeNbrs(X[i].size());
		for (size_t j=0, sz=X[i].size(); j<sz; ++j) W[i].SetNbr(j, X[i][j]);
	}
}

GalElement* Gda::GetGalElement(GeoDaWeight* w)
{
    GalElement *gal = 0;
    if (w->weight_type == GeoDaWeight::gal_type) {
        GalWeight *gal_w = dynamic_cast<GalWeight *>(w);
        gal = gal_w->gal;
    } else {
        GwtWeight *gwt_w = dynamic_cast<GwtWeight *>(w);
        GwtElement *gwt = gwt_w->gwt;
        gal = Gda::Gwt2Gal(gwt, gwt_w->num_obs);
    }
    return gal;
}

GalElement* Gda::NeighborMapToGal(std::vector<std::set<int> >& nbr_map)
{
	if (nbr_map.size() == 0) return 0;
	GalElement* gal = new GalElement[nbr_map.size()];
	if (!gal) return 0;
	for (int i=0, iend=nbr_map.size(); i<iend; i++) {
		gal[i].SetSizeNbrs(nbr_map[i].size());
		long cnt = 0;
		for (std::set<int>::iterator it=nbr_map[i].begin();
			 it != nbr_map[i].end(); it++) {
			gal[i].SetNbr(cnt++, *it);
		}
	}
	return gal;
}
