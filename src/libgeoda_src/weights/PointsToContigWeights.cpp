#include <limits>
#include <algorithm>
#include <map>
#include <string>
#include <sstream>
#include <set>
#include <list>

#include "PointsToContigWeights.h"

#define JC_VORONOI_IMPLEMENTATION
// If you wish to use doubles
//#define JCV_REAL_TYPE double
//#define JCV_ATAN2 atan2
//#define JCV_FLT_MAX 1.7976931348623157E+30

#include "../shape/jc_voronoi.h"

std::string jcv_point_str(const jcv_point& pt) 
{
    std::stringstream ss;
    ss << pt.x << "," << pt.y;
    return ss.str();
}

std::string jcv_edge_str(const jcv_graphedge* edge) 
{
    std::stringstream ss;
    jcv_real x0 = std::min(edge->pos[0].x, edge->pos[0].x);
    jcv_real x1 = std::max(edge->pos[1].x, edge->pos[1].x);
    jcv_real y0 = std::min(edge->pos[0].y, edge->pos[0].y);
    jcv_real y1 = std::max(edge->pos[1].y, edge->pos[1].y);
    ss << "[" << x0 << "," << y0 << "," << x1 << "," << y1 << "]";
    return ss.str();
}

/** If false returned, then an unexpected error.  Otherwise, neighbor map
 created successfully.  The presence of duplicates is indicated in
 duplicates_exists and the list of duplicates is filled in.
 */
bool gda::PointsToContiguity(const std::vector<double>& x,
    const std::vector<double>& y,
    bool queen,
    std::vector<std::set<int> >& nbr_map)
{
    size_t num_obs = x.size();

    double x_orig_min = std::numeric_limits<double>::max();
    double y_orig_min = std::numeric_limits<double>::max();
    double x_orig_max = std::numeric_limits<double>::min();
    double y_orig_max = std::numeric_limits<double>::min();

    for (size_t i=0; i<num_obs; ++i) {
        if (x_orig_min < x[i]) x_orig_min = x[i];
        if (x_orig_max > x[i]) x_orig_max = x[i];
        if (y_orig_min < y[i]) y_orig_min = y[i];
        if (y_orig_max > y[i]) y_orig_max = y[i];
    }

    double x_range = x_orig_max-x_orig_min;
    double y_range = y_orig_max-y_orig_min;

    // Add 2% offset to the bounding rectangle
	const double bb_pad = 0.02;
	// note data has been translated to origin and scaled
	double bb_xmin = x_orig_min - x_range * bb_pad;
	double bb_xmax = x_orig_max + x_range * bb_pad;
	double bb_ymin = y_orig_min - y_range * bb_pad;
	double bb_ymax = y_orig_max + y_range * bb_pad;

    // seed sites
    jcv_point* points = new jcv_point[num_obs];     
    for (size_t i=0; i< num_obs; ++i) {
        points[i].x = (float)x[i];
        points[i].y = (float)y[i];
    }

	nbr_map.clear();
	nbr_map.resize(num_obs);

    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));

    jcv_rect bounding_box;
    bounding_box.min.x = bb_xmin;
    bounding_box.min.y = bb_ymin;
    bounding_box.max.x = bb_xmax;
    bounding_box.max.y = bb_ymax;

    // create a voronoi graph
    jcv_diagram_generate(num_obs, (const jcv_point *)points, &bounding_box, &diagram);

    std::map<std::string, std::set<int> > edge_to_site;
    std::map<std::string, std::set<int> > jcvpoint_to_site;

    // edges
    const jcv_site* sites = jcv_diagram_get_sites(&diagram);

    for (size_t i=0; i<diagram.numsites; i++) {
        const jcv_site* site = &sites[i];
        const jcv_graphedge* e = sites[i].edges;

        while (e) {
            // shared edges will be visited by neighbor sites
            if (queen) {
                jcvpoint_to_site[jcv_point_str(e->pos[0])].insert(site->index);
                jcvpoint_to_site[jcv_point_str(e->pos[1])].insert(site->index);
            } else {
                edge_to_site[jcv_edge_str(e)].insert(site->index);
            }

            e = e->next;
        }
    }

    std::map<std::string, std::set<int> >::iterator it;
    std::set<int>::iterator nbr_it;
    if (queen) {
        for (it=jcvpoint_to_site.begin(); 
            it != jcvpoint_to_site.end(); ++it) 
        {
            // it->second are neighbors
            const std::set<int>& nbrs = it->second;
            for (nbr_it=nbrs.begin(); nbr_it!=nbrs.end(); ++nbr_it) {
                int ii = *nbr_it;
                for (nbr_it=nbrs.begin(); nbr_it!=nbrs.end(); ++nbr_it) {
                    int jj = *nbr_it;
                    if (ii != jj) {
                        nbr_map[ii].insert(jj);
                    }
                }
            }
        }
    } else {
        for (it=edge_to_site.begin(); 
            it != edge_to_site.end(); ++it) 
        {
            // it->second are neighbors
            const std::set<int>& nbrs = it->second;
            for (nbr_it=nbrs.begin(); nbr_it!=nbrs.end(); ++nbr_it) {
                int ii = *nbr_it;
                for (nbr_it=nbrs.begin(); nbr_it!=nbrs.end(); ++nbr_it) {
                    int jj = *nbr_it;
                    if (ii != jj) {
                        nbr_map[ii].insert(jj);
                    }
                }
            }
        }
    }
    delete[] points;
    jcv_diagram_free( &diagram );
    return true;
}

