#ifndef __JSGEODA_GAL_WEIGHT_H__
#define __JSGEODA_GAL_WEIGHT_H__

#include <vector>
#include <map>
#include <set>
#include "GeodaWeight.h"

class GalElement {
public:
	GalElement();
	void SetSizeNbrs(size_t sz, bool is_gal=false);
	void SetNbr(size_t pos, long n);
	void SetNbr(size_t pos, long n, double w);
	void SetNbrs(const GalElement& gal);
	const std::vector<long>& GetNbrs() const;
	const std::vector<double>& GetNbrWeights() const;
	void SortNbrs();
    void ReverseNbrs();
	long Size() const { return nbr.size(); }
	long operator[](size_t n) const { return nbr[n]; }
	double SpatialLag(const std::vector<double>& x) const;
	double SpatialLag(const double* x) const;
	double SpatialLag(const std::vector<double>& x, const int* perm) const;
    double GetRW(int idx);
    bool   Check(long nbrIdx);
   
    bool is_nbrAvgW_empty;
    std::vector<double> nbrAvgW;
    std::map<long, int> nbrLookup; // nbr_id, idx_in_nbrWeight
    
    void Update(const std::vector<bool>& undefs);

    int idx;

private:
	std::vector<long> nbr;
	std::vector<double> nbrWeight;
};

class GalWeight : public GeoDaWeight {
public:
	GalElement* gal;
    
	GalWeight() : gal(0) { weight_type = gal_type; }
    
	GalWeight(const GalWeight& gw);
    
	virtual ~GalWeight() { if (gal) delete [] gal; gal = 0; }
    
	static bool HasIsolates(GalElement *gal, int num_obs);
    
	virtual GalWeight& operator=(const GalWeight& gw);
    
	virtual bool HasIsolations() { return HasIsolates(gal, num_obs); }
    
    virtual void Update(const std::vector<bool>& undefs);
    
    virtual bool CheckNeighbor(int obs_idx, int nbr_idx);

    virtual const std::vector<long> GetNeighbors(int obs_idx);

    virtual const  std::vector<double> GetNeighborWeights(int obs_idx);

    virtual void GetNbrStats();
    
    virtual int GetNbrSize(int obs_idx);
    
    virtual double SpatialLag(int obs_idx,
                              const std::vector<double>& data);
    virtual bool   Save(const char* ofname,
                              const char* layer_name,
                              const char* id_var_name,
                              const std::vector<int>& id_vec);

    virtual bool   Save(const char* ofname,
                              const char* layer_name,
                              const char* id_var_name,
                              const std::vector<const char*>& id_vec);
};

namespace Gda {
	void MakeHigherOrdContiguity(size_t distance, size_t obs, GalElement* W, bool cummulative);
    GalElement* GetGalElement(GeoDaWeight* w);
    GalElement* NeighborMapToGal(std::vector<std::set<int> >& nbr_map);
}

#endif
