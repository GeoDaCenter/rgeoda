#include "GeodaWeight.h"

GeoDaWeight::GeoDaWeight(const GeoDaWeight& gw)
{
	GeoDaWeight::operator=(gw);
}

const GeoDaWeight& GeoDaWeight::operator=(const GeoDaWeight& gw)
{
	weight_type = gw.weight_type;
	wflnm = gw.wflnm;
	id_field = gw.id_field;
	title = gw.title;
	symmetry_checked = gw.symmetry_checked;
	is_symmetric = gw.is_symmetric;
	num_obs = gw.num_obs;
	sparsity = gw.sparsity;
	density = gw.density;
	min_nbrs = gw.min_nbrs;
	max_nbrs = gw.max_nbrs;
	mean_nbrs = gw.mean_nbrs;
	median_nbrs = gw.median_nbrs;
	is_internal_use = gw.is_internal_use;
	uid = gw.uid;
	
	return *this;
}

std::string GeoDaWeight::GetTitle()  const
{
	return title;
}

bool GeoDaWeight::IsSymmetric() const
{
    return is_symmetric;
}
double GeoDaWeight::GetSparsity() const
{
    return sparsity;
}
double GeoDaWeight::GetDensity() const
{
    return density;
}
int GeoDaWeight::GetNumObs() const
{
    return num_obs;
}

int GeoDaWeight::GetMinNbrs() const
{
    return min_nbrs;
}
int GeoDaWeight::GetMaxNbrs() const
{
    return max_nbrs;
}
double GeoDaWeight::GetMeanNbrs() const
{
    return mean_nbrs;
}
double GeoDaWeight::GetMedianNbrs() const
{
    return median_nbrs;
}
