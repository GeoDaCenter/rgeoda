#ifndef __JSGEODA_SPATIAL_IND_ALGS_H__
#define __JSGEODA_SPATIAL_IND_ALGS_H__

#include <list>
#include <set>
#include <sstream>
#include <vector>

#include "SpatialIndTypes.h"
#include "weights/GwtWeight.h"

#ifndef M_PI
#define M_PI            3.14159265358979310
#endif

namespace SpatialIndAlgs {
    
void to_3d_centroids(const std::vector<pt_2d>& pt2d,
										 std::vector<pt_3d>& pt3d);
void to_3d_centroids(const std::vector<pt_lonlat>& ptll,
										 std::vector<pt_3d>& pt3d);
void default_test();
void print_rtree_stats(rtree_box_2d_t& rtree);
void query_all_boxes(rtree_box_2d_t& rtree);
void knn_query(const rtree_pt_2d_t& rtree, int nn=6);
/** Will call more specialized knn_build as needed.  This routine will
 build the correct type of rtree automatically.  If is_arc false,
 then Euclidean distance is used and x, y are normal coordinates and
 is_mi ignored.  If is_arc is true, then arc distances are used and distances
 reported in either kms or miles according to is_mi. */

void apply_kernel(const GwtWeight* Wp, const std::string& kernel, bool use_kernel_diagnals = false);
GwtWeight* knn_build(const std::vector<gda::PointContents*>& points,
                     int nn,
                     bool is_arc, bool is_mi,
                     bool is_inverse=false, double power=1,
                     const std::string& kernel = "",
                     double bandwidth = 0,
                     bool adaptive_bandwidth = false,
                     bool use_kernel_diagnals = false);
GwtWeight* knn_build(const rtree_pt_2d_t& rtree, int nn=6,
                     bool is_inverse=false, double power=1,
                     const std::string& kernel = "",
                     double bandwidth = 0,
                     bool adaptive_bandwidth = false,
                     bool use_kernel_diagnals = false);
GwtWeight* knn_build(const rtree_pt_3d_t& rtree, int nn=6,
					 bool is_arc=false, bool is_mi=true,
                     bool is_inverse=false, double power=1,
                     const std::string& kernel = "",
                     double bandwidth = 0,
                     bool adaptive_bandwidth = false,
                     bool use_kernel_diagnals = false);
double est_thresh_for_num_pairs(const rtree_pt_2d_t& rtree, double num_pairs);
double est_thresh_for_avg_num_neigh(const rtree_pt_2d_t& rtree, double avg_n);
double est_avg_num_neigh_thresh(const rtree_pt_2d_t& rtree, double th,
								size_t trials=100);
/** If is_arc true, result is returned as radians. If
 * (x.size()*x.size()-1)/2 <= max_iters, exact value is returned */
double est_mean_distance(const std::vector<double>& x,
						 const std::vector<double>& y,
						 bool is_arc, size_t max_iters=300000);
/** If is_arc true, result is returned as radians. If
 * (x.size()*x.size()-1)/2 <= max_iters, exact value is returned*/
double est_median_distance(const std::vector<double>& x,
						   const std::vector<double>& y,
						   bool is_arc, size_t max_iters=300000);
/** Will call more specialized thresh_build as needed.  This routine will
 build the correct type of rtree automatically.  If is_arc false,
 then Euclidean distance is used and x, y are normal coordinates and
 is_mi ignored.  If is_arc is true, then arc distances are used and distances
 reported in either kms or miles according to is_mi.  When is_arc is true,
 the threshold input parameter is assumed to be in earth arc miles or kms
 according to is_mi. */
GwtWeight* thresh_build(const std::vector<double>& x,
                        const std::vector<double>& y,
                        double th,
                        double power,
                        bool is_arc, bool is_mi,
                        const std::string& kernel="", bool use_kernel_diagnals=false);
GwtWeight* thresh_build(const rtree_pt_2d_t& rtree, double th, double power
                        , const std::string& kernel="", bool use_kernel_diagnals=false);
double est_avg_num_neigh_thresh(const rtree_pt_3d_t& rtree, double th,
								size_t trials=100);
/** threshold th is the radius of intersection sphere with
  respect to the unit shpere of the 3d point rtree */
GwtWeight* thresh_build(const rtree_pt_3d_t& rtree, double th, double power,
                        bool is_mi, const std::string& kernel="",
                        bool use_kernel_diagnals=false);
/** Find the nearest neighbor for all points and return the maximum
 distance of all of these nearest neighbor pairs.  This is the minimum
 threshold distance such that all points have at least one neighbor.
 is_mi only relevant when is_arc is true.*/
double find_max_1nn_dist(const std::vector<double>& x,
						const std::vector<double>& y,
						bool is_arc, bool is_mi);
void get_pt_rtree_stats(const rtree_pt_2d_t& rtree,
						double& min_d_1nn, double& max_d_1nn,
						double& mean_d_1nn, double& median_d_1nn);
/** results returned in radians */
void get_pt_rtree_stats(const rtree_pt_3d_t& rtree,
						double& min_d_1nn, double& max_d_1nn,
						double& mean_d_1nn, double& median_d_1nn);
GwtWeight* knn_build(const rtree_pt_lonlat_t& rtree, int nn=6);
bool write_gwt(const GwtWeight* W, const std::string& layer_name,
			   const std::string& ofname, const std::string& vname,
			   const std::vector<int>& id_vec);
void fill_pt_rtree(rtree_pt_2d_t& rtree,
				   const std::vector<pt_2d>& pts);
void fill_pt_rtree(rtree_pt_lonlat_t& rtree,
				   const std::vector<pt_lonlat>& pts);
void fill_pt_rtree(rtree_pt_3d_t& rtree,
				   const std::vector<pt_3d>& pts);
struct LonLatPt {
	LonLatPt() : lon(0), lat(0) {}
	LonLatPt(double lon_, double lat_) : lon(lon_), lat(lat_) {}
	double lon;
	double lat;
};
std::ostream& operator<< (std::ostream &out, const LonLatPt& pt);
struct XyzPt {
	XyzPt() : x(0), y(0), z(0) {}
	XyzPt(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
	double x;
	double y;
	double z;
};
std::ostream& operator<< (std::ostream &out, const XyzPt& pt);

}
	
#endif

