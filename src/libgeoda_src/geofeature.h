#ifndef __GEODA_CENTER_SHP_FILE_H__
#define __GEODA_CENTER_SHP_FILE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <float.h>

namespace gda {
	
	enum ShapeType {
		NULL_SHAPE = 0, POINT_TYP = 1, POLY_LINE = 3, POLYGON = 5,
		MULTI_POINT = 8, POINT_Z = 11, POLY_LINE_Z = 13, POLYGON_Z = 15,
		MULTI_POINT_Z = 18, POINT_M = 21, POLY_LINE_M = 23, POLYGON_M = 25,
		MULTI_POINT_M = 28, MULTI_PATCH = 31
	};
	
	struct Point {
		Point() : x(0), y(0) {}
		Point(double x_s, double y_s) : x(x_s), y(y_s) {}
		double x;
		double y;
		bool equals(const Point* p, double precision_threshold=0.0) const {
			return (abs(x - p->x) <= precision_threshold &&
							abs(y - p->y) <= precision_threshold);
		}
		bool equals(const Point& p, double precision_threshold=0.0) const {
			return (abs(x - p.x) <= precision_threshold &&
							abs(y - p.y) <= precision_threshold);
		}
		double distance(const Point& p) const {
			return sqrt( (x - p.x)*(x-p.x) + (y-p.y)*(y-p.y) );
		}
	};
	
	bool operator==(Point const& a, Point const& b);
	
	struct LineSegment {
		LineSegment() : a(0,0), b(1,1) {}
		LineSegment(const Point& a_s, const Point& b_s, bool ordered=false)
		: a(a_s), b(b_s)
		{
			if (!ordered && ((a_s.x > b_s.x) || (a_s.x == b_s.x && a_s.y > b_s.y ))) {
				a = b_s; b = a_s;
			}
		}
		Point a;
		Point b;
	};
	
	struct GeometryContent {
		GeometryContent() : shape_type(0) {}
		virtual ~GeometryContent() {}
		GeometryContent(gda::ShapeType st) : shape_type(st) {}
		int shape_type; // byte 0, LE, one value from ShapeType enum
	};
	
	struct NullShapeContents : public GeometryContent {
		NullShapeContents() : GeometryContent(gda::NULL_SHAPE) {}
		virtual ~NullShapeContents() {}
		//shape_type = 0, byte 0, LE
	};
	
	struct PointContents : public GeometryContent {
		PointContents() : GeometryContent(gda::POINT_TYP), x(0), y(0) {}
		virtual ~PointContents() {}
		//shape_type = 1, byte 0, LE
		double x; // byte 4, LE
		double y; // byte 12, LE
	};
	
	struct PolyLineContents : public GeometryContent {
		PolyLineContents() : GeometryContent(gda::POLY_LINE), box(4,0),
		num_parts(0), num_points(0), parts(0), points(0) {}
		virtual ~PolyLineContents() {}
		//shape_type = 3, byte 0, LE
		std::vector<double> box; // byte 4, LE
		int num_parts; // byte 36, LE
		int num_points; // byte 40, LE
		std::vector<int> parts; // byte 44, array of size num_parts, LE
		// stores the first index in array for each part
		//Point points; // byte 44 + 4*num_parts, array of size num_parts, LE
		std::vector<Point> points;
	};
	
	struct PolygonContents : public GeometryContent {
		PolygonContents() : GeometryContent(gda::POLYGON), box(4,0),
		num_parts(0), num_points(0), parts(0), points(0) {}
		virtual ~PolygonContents() {}
		//shape_type = 5, byte 0, LE
		std::vector<double> box; // byte 4, LE minx, miny, maxx, maxy
		int num_parts; // byte 36, LE, the number of rings in the polygon
		int num_points; // byte 40, LE, total number of points for all rings
		std::vector<int> parts; // byte 44, array of size num_parts, LE
		std::vector<bool> holes;
		// stores the first index in array for each part
		std::vector<Point> points; // byte 44 + 4*num_parts, array of
		// size num_parts, LE
		
		bool intersect(PolygonContents* shp) { 
			return !(shp->box[0] > box[2] || shp->box[1] > box[3] ||
					 shp->box[2] < box[0] || shp->box[3] < box[1]);
		}
	};

	struct MainMap {
		MainMap() : num_obs(0), shape_type(NULL_SHAPE),
		bbox_x_min(DBL_MAX), bbox_y_min(DBL_MAX), bbox_x_max(DBL_MIN), bbox_y_max(DBL_MIN) {}
		virtual ~MainMap() {
			for (size_t i=0; i<records.size(); i++) {
				if(records[i]) delete records[i];
			}
			records.clear();
		}
		void set_bbox(double x, double y) {
			if ( x < bbox_x_min ) bbox_x_min = x;
            if ( x >= bbox_x_max ) bbox_x_max = x;
            if ( y < bbox_y_min ) bbox_y_min = y;
            if ( y >= bbox_y_max ) bbox_y_max = y;
		}

		int num_obs;
		gda::ShapeType shape_type;
		double bbox_x_min; // bounding box for X. similar below
		double bbox_y_min; // 
		double bbox_x_max; // 
		double bbox_y_max; // 
		std::vector<GeometryContent*> records;
	};
}

#endif
