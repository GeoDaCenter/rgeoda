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

#ifndef __GEODA_CENTER_GEODA_WEIGHTS_H__
#define __GEODA_CENTER_GEODA_WEIGHTS_H__

#include <vector>
#include <wx/string.h>

class GeoDaWeight {
public:
	GeoDaWeight() : symmetry_checked(false), num_obs(0) {}
	GeoDaWeight(const GeoDaWeight& gw);
    
	virtual ~GeoDaWeight() {}
	
    // following functions implemented in inherited classes:
    // GalWeights and GwtWeights
    virtual bool SaveDIDWeights(int num_obs,
                                std::vector<wxInt64>& newids,
                                std::vector<wxInt64>& stack_ids,
                                const wxString& ofname) = 0;
    
    virtual bool SaveSpaceTimeWeights(const wxString& ofname,
                                      const std::vector<wxString>& id_vec,
                                      const std::vector<wxString>& time_ids) = 0;

    virtual bool CheckNeighbor(int obs_idx, int nbr_idx) = 0;

    virtual const std::vector<long> GetNeighbors(int obs_idx) = 0;

    virtual void   Update(const std::vector<bool>& undefs) = 0;
    virtual bool   HasIsolates() = 0;
    virtual void   GetNbrStats() = 0;
    
    virtual int    GetNbrSize(int obs_idx) = 0;
    virtual double SpatialLag(int obs_idx,
                              const std::vector<double>& data) = 0;
    
    virtual bool   SaveToFile(const wxString& ofname,
                              const wxString& layer_name,
                              const wxString& id_var_name,
                              const std::vector<wxInt64>& id_vec) = 0;
    
    virtual bool   SaveToFile(const wxString& ofname,
                              const wxString& layer_name,
                              const wxString& id_var_name,
                              const std::vector<wxString>& id_vec) = 0;
    
    // functions:
    virtual double GetSparsity() const;
    virtual double GetDensity() const;
    virtual int    GetMinNumNbrs() const;
    virtual int    GetMaxNumNbrs() const;
    virtual double GetMeanNumNbrs() const;
    virtual double GetMedianNumNbrs() const;
    virtual int    GetNumObs() const;
    virtual bool   IsInternalUse() const { return is_internal_use; }
    
    // Others
    virtual const GeoDaWeight& operator=(const GeoDaWeight& gw);
   
    virtual wxString GetTitle() const; // returns portion of wflnm if title empty
   
    virtual wxString GetIDName() const { return id_field;}

    // Properties
	enum WeightType { gal_type, gwt_type };
	WeightType weight_type;
	wxString   wflnm; // filename
    wxString   id_field;
	wxString   title; // optional title.  Use wflnm if empty
	bool       symmetry_checked; // indicates validity of is_symmetric bool
	bool       is_symmetric; // true iff matrix is symmetric
	int        num_obs;
    double     sparsity;
    double     density;
    int        min_nbrs;
    int        max_nbrs;
    double     mean_nbrs;
    double     median_nbrs;
    bool       is_internal_use; // if internally used weights structure, will not be shown and used by users
};

#endif

