#ifndef __JSGEODA_CENTROID__
#define __JSGEODA_CENTROID__


// Assemble and include the file manually:
// 1. 12.0\VC\bin\x86_amd64\ml64.exe" /c ttmathuint_x86_64_msvc.asm (inside the ttmath folder).
// 2. Include the object file in compilation: ttmath/ttmathuint_x86_64_msvc.obj
#ifdef _WIN32
// or Disable assembly by adding the following line before #include ttmath.h
#define TTMATH_NOASM 1
#endif

#include <cmath>
#include "ttmath/ttmath.h"
#include "../geofeature.h"

/// Usage: `ttmath::Big<exponent, mantissa>`
typedef ttmath::Big<TTMATH_BITS(32), TTMATH_BITS(128)> DD;
//typedef ttmath::Big<TTMATH_BITS(64), TTMATH_BITS(128)> DD;
//typedef ttmath::Big<TTMATH_BITS(32), TTMATH_BITS(256)> DD;
//typedef ttmath::Big<TTMATH_BITS(64), TTMATH_BITS(256)> DD;

double DP_SAFE_EPSILON =  1e-15;

enum {
    RIGHT = -1,
    LEFT = 1,
    STRAIGHT = 0,
    FAILURE = 2
};

inline int
OrientationDD(DD const& dd)
{
    static DD const zero(0.0);
    if(dd < zero) {
        return RIGHT;
    }

    if(dd > zero) {
        return LEFT;
    }

    return STRAIGHT;
}

static int
orientation(double x)
{
    if(x < 0) {
        return RIGHT;
    }
    if(x > 0) {
        return LEFT;
    }
    return STRAIGHT;
}

class Orientation {
    public:

    /* A value that indicates an orientation or turn */
    enum {
        CLOCKWISE = -1,
        COLLINEAR = 0,
        COUNTERCLOCKWISE = 1,
        RIGHT = -1,
        LEFT = 1,
        STRAIGHT = 0
    };

    static int orientationIndexFilter(const gda::Point& pa, const gda::Point& pb, const gda::Point& pc)
    {
        double detsum;
        double const detleft = (pa.x - pc.x) * (pb.y - pc.y);
        double const detright = (pa.y - pc.y) * (pb.x - pc.x);
        double const det = detleft - detright;

        if(detleft > 0.0) {
            if(detright <= 0.0) {
                return orientation(det);
            }
            else {
                detsum = detleft + detright;
            }
        }
        else if(detleft < 0.0) {
            if(detright >= 0.0) {
                return orientation(det);
            }
            else {
                detsum = -detleft - detright;
            }
        }
        else {
            return orientation(det);
        }

        double const errbound = DP_SAFE_EPSILON * detsum;
        if((det >= errbound) || (-det >= errbound)) {
            return orientation(det);
        }
        return FAILURE;
    }

    static int index(const gda::Point& p1, const gda::Point& p2, const gda::Point& q)
    {
        // fast filter for orientation index
        // avoids use of slow extended-precision arithmetic in many cases
        int index = orientationIndexFilter(p1, p2, q);
        if(index <= 1) {
            return index;
        }

        // normalize coordinates
        DD dx1 = DD(p2.x) + DD(-p1.x);
        DD dy1 = DD(p2.y) + DD(-p1.y);
        DD dx2 = DD(q.x) + DD(-p2.x);
        DD dy2 = DD(q.y) + DD(-p2.y);

        // sign of determinant - inlined for performance
        DD mx1y2(dx1 * dy2);
        DD my1x2(dy1 * dx2);
        DD d = mx1y2 - my1x2;
        return OrientationDD(d);
    }

    // detect orientation of a ring
    static bool isCCW(const std::vector<gda::Point>& pts, int start, int end)
    {
        // # of points without closing endpoint
        int nPts = end - start;

        if (nPts < 3) return false;

        // find highest point
        const gda::Point* hiPt = &pts[start];
        int hiIndex = start;

        for(int i = start+1; i <= end; ++i) {
            const gda::Point* p = &pts[i];
            if(p->y > hiPt->y) {
                hiPt = p;
                hiIndex = i;
            }
        }
        // find distinct point before highest point
        int iPrev = hiIndex;
        do {
            if(iPrev == start) {
                iPrev = end;
            }
            iPrev = iPrev - 1;
        }
        while(pts[iPrev].equals(hiPt) && iPrev != hiIndex);

        // find distinct point after highest point
        size_t iNext = hiIndex;
        do {
            iNext = (iNext + 1 - start) % nPts;
            iNext = iNext + start;
        }
        while(pts[iNext].equals(hiPt) && iNext != hiIndex);

        const gda::Point* prev = &pts[iPrev];
        const gda::Point* next = &pts[iNext];

        /*
        * This check catches cases where the ring contains an A-B-A
        * configuration of points.
        * This can happen if the ring does not contain 3 distinct points
        * (including the case where the input array has fewer than 4 elements),
        * or it contains coincident line segments.
        */

        if(prev->equals(hiPt) || next->equals(hiPt) || prev->equals(next)) {
            return false;
            // MD - don't bother throwing exception,
            // since this isn't a complete check for ring validity
            //throw  IllegalArgumentException("degenerate ring (does not contain 3 distinct points)");
        }

        int disc = Orientation::index(*prev, *hiPt, *next);

        /**
         *  If disc is exactly 0, lines are collinear.
         * There are two possible cases:
         *  (1) the lines lie along the x axis in opposite directions
         *  (2) the lines lie on top of one another
         *
         *  (1) is handled by checking if next is left of prev ==> CCW
         *  (2) should never happen, so we're going to ignore it!
         *  (Might want to assert this)
         */
        bool isCCW = false;

        if(disc == 0) {
            // poly is CCW if prev x is right of next x
            isCCW = (prev->x > next->x);
        }
        else {
            // if area is positive, points are ordered CCW
            isCCW = (disc > 0);
        }

        return isCCW;
    }
};



class Centroid {

public:
    Centroid(gda::PolygonContents* poly)
        :
        areasum2(0.0),
        totalLength(0.0),
        ptCount(0)
    {
        // The first element in the array represents the exterior ring.
        // Any subsequent elements represent interior rings
        for (size_t p = 0; p < poly->num_parts; ++p) {
            int start = poly->parts[p];
            int end = p+1 < poly->num_parts ? poly->parts[p+1] : poly->num_points;
            if (poly->holes[p]) {
                addHole(poly, start, end-1);
            } else {
                addShell(poly, start, end -1);
            }
        }
    }

    bool getCentroid(gda::PointContents& cent) const
    {
        if(std::abs(areasum2) > 0.0) {
            cent.x = cg3.x / 3 / areasum2;
            cent.y = cg3.y / 3 / areasum2;
        }
        else if(totalLength > 0.0) {
            // if polygon was degenerate, compute linear centroid instead
            cent.x = lineCentSum.x / totalLength;
            cent.y = lineCentSum.y / totalLength;
        }
        else if(ptCount > 0) {
            cent.x = ptCentSum.x / ptCount;
            cent.y = ptCentSum.y / ptCount;
        }
        else {
            return false;
        }
        return true;
    }

private:

    gda::Point areaBasePt;
    gda::Point triangleCent3;
    gda::Point lineCentSum;
    gda::Point ptCentSum;
    gda::Point cg3;

    double areasum2;
    double totalLength;
    int ptCount;


    void setAreaBasePoint(gda::Point& basePt)
    {
        areaBasePt.x = basePt.x;
        areaBasePt.y = basePt.y;
    }

    void centroid3(const gda::Point& p1, const gda::Point& p2, const gda::Point& p3, gda::Point& c)
    {
        c.x = p1.x + p2.x + p3.x;
        c.y = p1.y + p2.y + p3.y;
    }

    double area2(const gda::Point& p1, const gda::Point& p2, const gda::Point& p3)
    {
        return
            (p2.x - p1.x) * (p3.y - p1.y) -
            (p3.x - p1.x) * (p2.y - p1.y);
    }

    void addPoint(const gda::Point& pt)
    {
        ptCount += 1;
        ptCentSum.x += pt.x;
        ptCentSum.y += pt.y;
    }

    void addHole(const gda::PolygonContents* poly, int start, int end)
    {
        // NOTE: shapefiles like ESRI Shapefile record the shells and holes in the same orientation
        // if read by some libraries e.g. gdal the orientation of holes will be different than shells
        // so the following line should be changed according to different orientattion style of input
        // datasource
        bool isPositiveArea = !Orientation::isCCW(poly->points, start, end);
        for(size_t i = start, e = end; i < e; ++i) {
            addTriangle(areaBasePt, poly->points[i], poly->points[i + 1], isPositiveArea);
        }
        addLineSegments(poly->points, start, end);
    }

    void addTriangle(const gda::Point& p0, const gda::Point& p1, const gda::Point& p2, bool isPositiveArea)
    {
        double sign = (isPositiveArea) ? 1.0 : -1.0;
        centroid3(p0, p1, p2, triangleCent3);
        double a2 = area2(p0, p1, p2);
        cg3.x += sign * a2 * triangleCent3.x;
        cg3.y += sign * a2 * triangleCent3.y;
        areasum2 += sign * a2;
    }

    void addLineSegments(const std::vector<gda::Point>& pts, int start, int end)
    {
        size_t npts = end - start + 1;
        double lineLen = 0.0;
        for(size_t i = start; i < end; i++) {
            double segmentLen = pts[i].distance(pts[i + 1]);
            if(segmentLen == 0.0) {
                continue;
            }

            lineLen += segmentLen;

            double midx = (pts[i].x + pts[i + 1].x) / 2;
            lineCentSum.x += segmentLen * midx;
            double midy = (pts[i].y + pts[i + 1].y) / 2;
            lineCentSum.y += segmentLen * midy;
        }
        totalLength += lineLen;
        if(lineLen == 0.0 && npts > 0) {
            addPoint(pts[start]);
        }
    }

    void addShell(gda::PolygonContents* poly, int start, int end)
    {
        // ext ring
        size_t len = end - start + 1;
        if(len > 0) {
            setAreaBasePoint(poly->points[start]);
        }
        bool isPositiveArea = ! Orientation::isCCW(poly->points, start, end);
        for(size_t i = start; i < end; ++i) {
            addTriangle(areaBasePt, poly->points[i], poly->points[i + 1], isPositiveArea);
        }
        addLineSegments(poly->points, start, end);
    }

};
#endif
