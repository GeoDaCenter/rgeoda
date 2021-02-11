#include <iostream>
#include <fstream>
#include <sstream>

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

const char* WeightUtils::ReadIdField(const char* w_fname)
{
  /*
  std::ifstream file(w_fname);
  if (!(file.is_open() && file.good())) return "";

  // Header line is identical for GWT and GAL
  // First determine if header line is correct
  // Can be either: int int string string  (type n_obs filename field)
  // or : int (n_obs)
  std::string str;
  std::getline(file, str);
  std::stringstream ss(str, std::stringstream::in | std::stringstream::out);

  std::string line;
  std::getline(ss, line);
  std::string header(line);

  // detect if header contains string with empty space, which should be quoted
  if (header.find("\""))
   */
  return 0;
}
