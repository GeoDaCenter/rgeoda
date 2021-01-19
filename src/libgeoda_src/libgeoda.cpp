
#include <iostream>
#include <string.h>
#include <float.h>

#ifdef _MSC_VER
#include <cstddef>
#endif

#include "./shapelib/shapefil.h"
#include "./weights/GeodaWeight.h"
#include "./pg/geoms.h"
#include "./pg/utils.h"
#include "./shape/centroid.h"
#include "geofeature.h"
#include "gda_interface.h"

#include "libgeoda.h"

const std::string GeoDa::DT_STRING = "string";
const std::string GeoDa::DT_INTEGER= "integer";
const std::string GeoDa::DT_NUMERIC = "numeric";

GeoDa* CreateGeoDaFromGPD(const std::string& layer_name,
                         const std::string& map_type,
                         const std::vector<unsigned char> &wkbs,
                         const std::vector<int>& wkb_bytes_len)
{
    return new GeoDa(layer_name, map_type, wkbs, wkb_bytes_len);
}


GeoDa* CreateGeoDaFromSHP(const char* pDsPath, const char* layer_name)
{
    return new GeoDa(pDsPath, layer_name);
}

GeoDaColumn* ToGeoDaColumn(GeoDaStringColumn* col)
{
    return dynamic_cast<GeoDaColumn*>(col);
}

GeoDaColumn *ToGeoDaColumn(GeoDaIntColumn *col) {
    return dynamic_cast<GeoDaColumn*>(col);
}

GeoDaColumn *ToGeoDaColumn(GeoDaRealColumn *col) {
    return dynamic_cast<GeoDaColumn*>(col);
}

int test() { return 100;}

// this constructor is for R
GeoDa::GeoDa(const std::string &layer_name,
             const std::string& map_type,
             int num_features,
             unsigned char* wkbs,
             const std::vector<int>& wkb_bytes_len)
: numObs(num_features), numCols(0), table(NULL)
{
    main_map = new gda::MainMap();
    Init(layer_name, map_type, num_features, wkbs, wkb_bytes_len);
}

// this constructor is for Python
GeoDa::GeoDa(const std::string& layer_name,
             const std::string& map_type,
             const std::vector<unsigned char> &wkbs,
             const std::vector<int>& wkb_bytes_len)
: numObs(wkb_bytes_len.size()), numCols(0), table(NULL)
{
    main_map = new gda::MainMap();
	Init(layer_name, map_type, wkb_bytes_len.size(), (unsigned char*)(&wkbs[0]), wkb_bytes_len);
}

// this constructor is for reading from ESRI shapefile
GeoDa::GeoDa(const char* poDsPath, const char* layer_name)
: numLayers(0), numObs(0)
{
    main_map = new gda::MainMap();
    table = new GeoDaTable();

    char dbfpath[512];
    strcpy(dbfpath, poDsPath);
    strncpy(dbfpath + strlen(poDsPath)-3, "dbf", 4);

    ReadShapefile(poDsPath);
    ReadDbffile(dbfpath);
}

GeoDa::~GeoDa() {
    if (main_map) {
        delete main_map;
    }
    for (size_t i=0; i<centroids.size(); ++i) {
        if (centroids[i]) {
            delete centroids[i];
        }
    }
}

void GeoDa::Init(const std::string &layer_name,
                 const std::string& map_type,
                 int num_features,
                 unsigned char* wkbs,
                 const std::vector<int>& wkb_bytes_len)
{
    if (map_type == "map_polygons") {
        //mapType = polygon_type;
        main_map->shape_type = gda::POLYGON;
    } else if (map_type == "map_points") {
        main_map->shape_type = gda::POINT_TYP;
    } else if (map_type == "map_lines") {
        main_map->shape_type = gda::POLY_LINE;
    } else {
#ifndef __RGEODA__
        std::cout << "map type is not supported." << std::endl;
#endif
    }

    main_map->num_obs = num_features;

    // create features
    unsigned long long wkb_offset = 0;
    for (int i=0; i<num_features; ++i) {
        LWGEOM *lwgeom = lwgeom_from_wkb(wkbs + wkb_offset, wkb_bytes_len[i], LW_PARSER_CHECK_ALL);
        wkb_offset += wkb_bytes_len[i];

        if (lwgeom_is_empty(lwgeom)) {
            //lwdebug(4, "build_pg_geoad: addnullgeometry()");
            AddNullGeometry();
        } else {
            if (lwgeom->type == POINTTYPE) {
                //lwdebug(4, "build_pg_geoda: add point");
                LWPOINT *pt = lwgeom_as_lwpoint(lwgeom);
                AddPoint(pt);
            } else if (lwgeom->type == MULTIPOINTTYPE) {
                //lwdebug(4, "build_pg_geoda: add multi point");
                LWMPOINT *mpt = lwgeom_as_lwmpoint(lwgeom);
                AddMultiPoint(mpt);
            } else if (lwgeom->type == POLYGONTYPE) {
                //lwdebug(4, "build_pg_geoda: add polygon");
                LWPOLY *poly = lwgeom_as_lwpoly(lwgeom);
                AddPolygon(poly);
            } else if (lwgeom->type == MULTIPOLYGONTYPE) {
                //lwdebug(4, "build_pg_geoda: add multi-polygon");
                LWMPOLY *mpoly = lwgeom_as_lwmpoly(lwgeom);
                AddMultiPolygon(mpoly);
            } else {
                //lwdebug(4, "Unknown WKB type %s\n", lwtype_name(lwgeom->type));
                AddNullGeometry();
            }
        }
        /* Free both the original and geometries */
        lwgeom_free(lwgeom);
    }
}


void GeoDa::AddNullGeometry() {
    this->main_map->records.push_back(new gda::NullShapeContents());
}

void GeoDa::AddPoint(LWPOINT *lw_pt) {
    /* Grab the point: note getPoint4d will correctly handle
	the case where the POINTs don't contain Z or M coordinates */
    POINT4D p4d;
    p4d = getPoint4d(lw_pt->point, 0);

    gda::PointContents* pt = new gda::PointContents();
    pt->x = p4d.x;
    pt->y = p4d.y;
    this->main_map->set_bbox(pt->x,  pt->y);
    this->main_map->records.push_back(pt);
}

void GeoDa::AddMultiPoint(LWMPOINT *lw_mpt) {
    /* Grab the points: note getPoint4d will correctly handle
	the case where the POINTs don't contain Z or M coordinates */
    POINT4D p4d;
    for (size_t i = 0; i < lw_mpt->ngeoms; i++)
    {
        p4d = getPoint4d(lw_mpt->geoms[i]->point, 0);

        gda::PointContents* pt = new gda::PointContents();
        pt->x = p4d.x;
        pt->y = p4d.y;
        this->main_map->set_bbox(pt->x,  pt->y);
        this->main_map->records.push_back(pt);
        break;// only take the first point, even it has multipoints
    }
}

void GeoDa::AddPolygon(LWPOLY *lw_poly) {
    size_t  i, j;
    POINT4D p4d;
    /* Allocate storage for ring pointers */
    int shppointtotal = 0, shppoint = 0;

    /* First count through all the points in each ring so we now how much memory is required */
    for (i = 0; i < lw_poly->nrings; i++)
        shppointtotal += lw_poly->rings[i]->npoints;


    gda::PolygonContents *poly = new gda::PolygonContents();
    poly->num_parts = 0;
    poly->num_points = 0;

    double minx = DBL_MAX;
    double miny = DBL_MAX;
    double maxx = DBL_MIN;
    double maxy = DBL_MIN;
    double x, y;

    /* Iterate through each ring setting up shpparts to point to the beginning of each ring */
    for (i = 0; i < lw_poly->nrings; i++) {
        /* For each ring, store the integer coordinate offset for the start of each ring */
        poly->num_parts += 1;
        poly->parts.push_back(shppoint);

        bool is_hole = i > 0 ? true : false;
        poly->holes.push_back(is_hole);

        for (j = 0; j < lw_poly->rings[i]->npoints; j++) {
            p4d = getPoint4d(lw_poly->rings[i], j);

            x = p4d.x;
            y = p4d.y;

            poly->points.push_back(gda::Point(x,y));
            poly->num_points += 1;

            if ( x < minx ) minx = x;
            if ( x >= maxx ) maxx = x;
            if ( y < miny ) miny = y;
            if ( y >= maxy ) maxy = y;

            /* Increment the point counter */
            shppoint++;
        }

        /*
		 * First ring should be clockwise,
		 * other rings should be counter-clockwise
		 */
    }

    poly->box.resize(4);
    poly->box[0] = minx;
    poly->box[1] = miny;
    poly->box[2] = maxx;
    poly->box[3] = maxy;

    this->main_map->set_bbox(minx, miny);
    this->main_map->set_bbox(maxx, maxy);
    this->main_map->records.push_back(poly);
}

void GeoDa::AddMultiPolygon(LWMPOLY *lw_mpoly) {
    POINT4D p4d;
    uint32_t i, j, k;

    int shppointtotal = 0, shppoint = 0, shpringtotal = 0, shpring = 0;

    /* NOTE: Multipolygons are stored in shapefiles as Polygon* shapes with multiple outer rings */

    /* First count through each ring of each polygon so we now know much memory is required */
    for (i = 0; i < lw_mpoly->ngeoms; i++) {
        for (j = 0; j < lw_mpoly->geoms[i]->nrings; j++) {
            shpringtotal++;
            shppointtotal += lw_mpoly->geoms[i]->rings[j]->npoints;
        }
    }

    gda::PolygonContents *poly = new gda::PolygonContents();
    poly->num_parts = 0;
    poly->num_points = 0;

    double minx = DBL_MAX;
    double miny = DBL_MAX;
    double maxx = DBL_MIN;
    double maxy = DBL_MIN;
    double x, y;

    /* Iterate through each ring of each polygon in turn */
    for (i = 0; i < lw_mpoly->ngeoms; i++) {
        for (j = 0; j < lw_mpoly->geoms[i]->nrings; j++) {
            /* For each ring, store the integer coordinate offset for the start of each ring */
            poly->parts.push_back(shppoint);
            poly->num_parts += 1;

            bool is_hole = j > 0 ? true : false;
            poly->holes.push_back(is_hole);

            LWDEBUGF(4, "Ring offset: %d", shpring);

            for (k = 0; k < lw_mpoly->geoms[i]->rings[j]->npoints; k++) {
                p4d = getPoint4d(lw_mpoly->geoms[i]->rings[j], k);

                x = p4d.x;
                y = p4d.y;

                poly->points.push_back(gda::Point(x, y));
                poly->num_points += 1;

                if (x < minx) minx = x;
                if (x >= maxx) maxx = x;
                if (y < miny) miny = y;
                if (y >= maxy) maxy = y;

                /* Increment the point counter */
                shppoint++;
            }
            /*
			* First ring should be clockwise,
			* other rings should be counter-clockwise
			*/
        }
        /* Increment the ring counter */
        shpring++;
    }

    poly->box.resize(4);
    poly->box[0] = minx;
    poly->box[1] = miny;
    poly->box[2] = maxx;
    poly->box[3] = maxy;

    this->main_map->set_bbox(minx, miny);
    this->main_map->set_bbox(maxx, maxy);
    this->main_map->records.push_back(poly);
}


gda::MainMap& GeoDa::GetMainMap()
{
    return *this->main_map;
}

const std::vector<gda::PointContents*>& GeoDa::GetCentroids()
{
    // copy centroid from OGRGeometry
   if (centroids.empty()) {
       if (this->GetMapType() == gda::POINT_TYP) {
           this->centroids.resize(this->GetNumObs());
           for (size_t i=0; i<this->centroids.size(); ++i) {
               this->centroids[i] = new gda::PointContents;
               this->centroids[i]->x = ((gda::PointContents*)(this->main_map->records[i]))->x;
               this->centroids[i]->y = ((gda::PointContents*)(this->main_map->records[i]))->y;
           }
       } else if (this->GetMapType() == gda::POLYGON) {
           this->centroids.resize(this->GetNumObs());
           for (size_t i=0; i<this->centroids.size(); ++i) {
               gda::PolygonContents* poly = (gda::PolygonContents*)this->main_map->records[i];
               Centroid cent(poly);
               this->centroids[i] = new gda::PointContents;
               cent.getCentroid(*this->centroids[i]);
           }
       } else {
           lwerror("Enter PostGeoDa::GetCentroids() shape_type=%d not correct.", this->main_map->shape_type);
       }
    }
    return centroids;
}

std::vector<bool> GeoDa::GetNullValues(std::string col_name) {
    std::vector<bool> rst;
    /*
    if (fieldNameIdx.find(col_name) != fieldNameIdx.end()) {
        int iField = fieldNameIdx[col_name];
        OGRFeature *feature = NULL;
        for (size_t i=0; i<numObs; ++i) {
            feature = features[i];
            bool val = feature->IsFieldNull(iField);
            rst.push_back(val);
        }
    }
     */
    return rst;
}

std::vector<double> GeoDa::GetNumericCol(std::string col_name)
{
    std::vector<double> rst;
    if (table) {
        GeoDaColumn* col = table->GetColumn(col_name);
        if (col)  {
            if (col->field_type == GeoDaColumn::integer_type) {
                GeoDaIntColumn* r_col = dynamic_cast<GeoDaIntColumn*>(col);
                for (size_t i=0; i<r_col->data.size(); ++i) {
                    rst.push_back(r_col->data[i]);
                }

            } else if (col->field_type == GeoDaColumn::real_type) {
                GeoDaRealColumn* r_col = dynamic_cast<GeoDaRealColumn*>(col);
                rst = r_col->data;
            }
        }
    }
    return rst;
}

std::vector<long long> GeoDa::GetIntegerCol(std::string col_name) {
    std::vector<long long> rst;
    if (table) {
        GeoDaColumn* col = table->GetColumn(col_name);
        if (col)  {
            if (col->field_type == GeoDaColumn::integer_type) {
                GeoDaIntColumn* r_col = dynamic_cast<GeoDaIntColumn*>(col);
                rst = r_col->data;

            } else if (col->field_type == GeoDaColumn::real_type) {
                GeoDaRealColumn* r_col = dynamic_cast<GeoDaRealColumn*>(col);
                for (size_t i=0; i<r_col->data.size(); ++i) {
                    rst.push_back(r_col->data[i]);
                }
            }
        }
    }
    return rst;
}

std::vector<std::string> GeoDa::GetStringCol(std::string col_name) {
    std::vector<std::string> rst;
    if (table) {
        GeoDaColumn* col = table->GetColumn(col_name);
        if (col)  {
            if (col->field_type == GeoDaColumn::integer_type) {
                GeoDaIntColumn* r_col = dynamic_cast<GeoDaIntColumn*>(col);
                for (size_t i=0; i<r_col->data.size(); ++i) {
                    std::stringstream ss;
                    ss << r_col->data[i];
                    rst.push_back(ss.str());
                }

            } else if (col->field_type == GeoDaColumn::real_type) {
                GeoDaRealColumn* r_col = dynamic_cast<GeoDaRealColumn*>(col);
                for (size_t i=0; i<r_col->data.size(); ++i) {
                    std::stringstream ss;
                    ss << r_col->data[i];
                    rst.push_back(ss.str());
                }
            } else {
                GeoDaStringColumn* r_col = dynamic_cast<GeoDaStringColumn*>(col);
                rst = r_col->data;
            }
        }
    }
    return rst;
}

std::vector<std::string> GeoDa::GetFieldNames() {
    if (fieldNames.empty()) {
        if (table) {
            size_t n_cols = table->GetNumCols();
            for (size_t i=0; i<n_cols; ++i) {
                GeoDaColumn* col = table->GetColumn(i);
                fieldNames.push_back(col->name);
            }
        }
    }
    return fieldNames;
}

std::vector<std::string> GeoDa::GetFieldTypes() {
    if (fieldTypes.empty()) {
        if (table) {
            size_t n_cols = table->GetNumCols();
            for (size_t i=0; i<n_cols; ++i) {
                GeoDaColumn* col = table->GetColumn(i);
                if (col->field_type == GeoDaColumn::integer_type)
                    fieldTypes.push_back("integer");
                else if (col->field_type == GeoDaColumn::real_type)
                    fieldTypes.push_back("real");
                else
                    fieldTypes.push_back("string");
            }
        }
    }
    return fieldTypes;
}

int GeoDa::GetNumObs() const {
    return this->main_map->num_obs;
}

int GeoDa::GetNumCols() const {
    if (table) {
        return table->GetNumCols();
    }
    return 0;
}

int GeoDa::GetMapType()
{
    return main_map->shape_type;
}

void GeoDa::ReadDbffile(const char* fpath)
{
    DBFHandle hDBF = DBFOpen(fpath, "rb");
    if (hDBF == NULL) {
        // unable to open throw error
        return;
    }
    int		nWidth, nDecimals;
    char	szTitle[12];
    int     nRecords = DBFGetRecordCount(hDBF);

    for( int i = 0; i < DBFGetFieldCount(hDBF); i++ ) {
        DBFFieldType eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );

        std::vector<bool> undefs(nRecords, false);
        std::vector<long long> val_i;
        std::vector<double> val_d;
        std::vector<std::string> val_s;

        for( int iRecord = 0; iRecord < nRecords; iRecord++ ) {
            if (DBFIsAttributeNULL( hDBF, iRecord, i ) ) {
                undefs[iRecord] = true;
                continue;
            }
            if ( eType == FTInteger) {
                int val = DBFReadIntegerAttribute(hDBF, iRecord, i);
                val_i.push_back(val);
            } else if ( eType == FTDouble) {
                double val = DBFReadDoubleAttribute( hDBF, iRecord, i );
                val_d.push_back(val);
            } else {
                // others as FTString
                const char *val = DBFReadStringAttribute(hDBF, iRecord, i);
                val_s.push_back(val);
            }
        }

        if (eType == FTInteger)
            table->AddIntColumn(szTitle, val_i, undefs);
        else if (eType == FTDouble)
            table->AddRealColumn(szTitle, val_d, undefs);
        else
            table->AddStringColumn(szTitle, val_s, undefs);
    }
    DBFClose(hDBF);
}

void GeoDa::ReadShapefile(const char* fpath)
{


    SHPHandle hSHP = SHPOpen(fpath, "rb" );
    if( hSHP == NULL ) {
        // unable to open, throw error
        return;
    }
    int	nShapeType, nEntities, i, j;
    double 	adfMinBound[4], adfMaxBound[4];
    // get bounds
    SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );
    main_map->bbox_x_min = adfMinBound[0];
    main_map->bbox_y_min = adfMinBound[1];
    main_map->bbox_x_max = adfMaxBound[0];
    main_map->bbox_y_max = adfMaxBound[1];
    main_map->num_obs = nEntities;

    // read all shapes
    for( i = 0; i < nEntities; i++ ) {
        j = 0;
        SHPObject *psShape = SHPReadObject( hSHP, i );
        if( psShape == NULL ) {
            // can't read shape
            this->main_map->records.push_back(new gda::NullShapeContents());
            continue;
        }
        // bounds
        switch ( nShapeType ) {
            case SHPT_POINT: {
                main_map->shape_type = gda::POINT_TYP;
                gda::PointContents* pc = new gda::PointContents();
                pc->x = *psShape->padfX;
                pc->y = *psShape->padfY;
                main_map->records.push_back(pc);
                break;
            }

            case SHPT_MULTIPOINT: {
                main_map->shape_type = gda::POINT_TYP;
                gda::PointContents* pc = new gda::PointContents();
                //for( j=0; j<psShape->nVertices; ++j ) {
                    pc->x = psShape->padfX[0];
                    pc->y = psShape->padfY[0];
                    main_map->records.push_back(pc);
                    // only first point is used
                //}
                break;
            }

            case SHPT_ARC: {
                main_map->shape_type = gda::POLY_LINE;
                gda::PolyLineContents* pc = new gda::PolyLineContents();
                pc->num_parts = psShape->nParts;
                pc->num_points = psShape->nVertices;
                pc->box[0] = psShape->dfXMin;
                pc->box[1] = psShape->dfYMin;
                pc->box[2] = psShape->dfXMax;
                pc->box[3] = psShape->dfYMax;
                double x,y;
                for (j = 0; j < psShape->nParts; ++j) {
                    int itEnd = (j + 1 < psShape->nParts) ? psShape->panPartStart[j + 1] : psShape->nVertices;
                    pc->parts.push_back(psShape->panPartStart[j]);

                    for (int k = psShape->panPartStart[j]; k < itEnd; ++k) {
                        // line append
                        x = psShape->padfX[k];
                        y = psShape->padfY[k];
                        pc->points.push_back(gda::Point(x,y));
                    }
                }
                main_map->records.push_back(pc);
                break;
            }

            case SHPT_POLYGON: {
                main_map->shape_type = gda::POLYGON;
                gda::PolygonContents* pc = new gda::PolygonContents();
                pc->num_parts = psShape->nParts;
                pc->num_points = psShape->nVertices;
                pc->box[0] = psShape->dfXMin;
                pc->box[1] = psShape->dfYMin;
                pc->box[2] = psShape->dfXMax;
                pc->box[3] = psShape->dfYMax;
                double x,y;
                for (j = 0; j < psShape->nParts; ++j) {
                    int itStart = psShape->panPartStart[j];
                    int itEnd = (j + 1 < psShape->nParts) ? psShape->panPartStart[j + 1] : psShape->nVertices;
                    pc->parts.push_back(itStart);
                    pc->holes.push_back(j>0);
                    for (int k = itStart; k < itEnd; ++k) {
                        // ring append
                        x = psShape->padfX[k];
                        y = psShape->padfY[k];
                        pc->points.push_back(gda::Point(x,y));
                    }
                }
                main_map->records.push_back(pc);
                break;
            }
            break;
        }
    }
    SHPClose( hSHP );
}

// The following function will be reserved for working with GDAL
void GeoDa::ReadAllFeatures()
{
    /*
    bool has_null_geometry = false;

    // get geometry envelope
    //main_map->bbox_x_min = std::numeric_limits<double>::max();
    //main_map->bbox_y_min = std::numeric_limits<double>::max();
    //main_map->bbox_x_max = std::numeric_limits<double>::lowest();
    //main_map->bbox_y_max = std::numeric_limits<double>::lowest();

    // resize geometry records
    main_map->num_obs = numObs;


    //read OGR geometry features
    int feature_counter =0;
    for (int row_idx=0; row_idx < numObs; row_idx++ ) {
        OGRFeature* feature = features[row_idx];
        OGRGeometry* geometry= feature->GetGeometryRef();
        if (geometry == 0) {
            has_null_geometry = true;
            this->main_map->records.push_back(new gda::NullShapeContents());
            continue;
        }
        OGRwkbGeometryType eType = wkbFlatten(geometry->getGeometryType());
        // sometime OGR can't return correct value from GetGeomType() call
        if (eType == wkbPoint) {
            gda::PointContents* pc = new gda::PointContents();
            OGRPoint* p = (OGRPoint *) geometry;
            pc->x = p->getX();
            pc->y = p->getY();
            main_map->set_bbox(pc->x, pc->y);
            main_map->records.push_back(pc);

        } else if (eType == wkbMultiPoint) {
            gda::PointContents* pc = new gda::PointContents();
            OGRMultiPoint* mp = (OGRMultiPoint*) geometry;
            int n_geom = mp->getNumGeometries();
            for (size_t i = 0; i < n_geom; i++ ) {
                // only consider first point
                OGRGeometry* ogrGeom = mp->getGeometryRef(i);
                OGRPoint* p = (OGRPoint*)ogrGeom;
                pc->x = p->getX();
                pc->y = p->getY();
                main_map->set_bbox(pc->x, pc->y);
                //if (noExtent) GetExtent(main_map, pc, row_idx);
            }
            main_map->records.push_back(pc);

        } else if (eType == wkbPolygon || eType == wkbCurvePolygon ) {
            gda::PolygonContents* pc = new gda::PolygonContents();

            OGRPolygon* p = (OGRPolygon *) geometry;
            CopyEnvelope(p, pc);

            OGRLinearRing* pLinearRing = 0;
            int numPoints= 0;

            // interior rings
            int ni_rings = p->getNumInteriorRings();
            // resize parts memory, 1 is for exterior ring,
            pc->num_parts = ni_rings + 1;
            for (size_t j=0; j < pc->num_parts; j++ ) {
                pLinearRing = j==0 ? p->getExteriorRing() : p->getInteriorRing(j-1);
                pc->holes.push_back(j>0);
                pc->parts.push_back(numPoints);
                numPoints += pLinearRing->getNumPoints();
            }
            // resize points memory
            pc->num_points = numPoints;
            // read points
            double x, y;
            for (size_t j=0; j < pc->num_parts; j++) {
                pLinearRing = j==0 ? p->getExteriorRing() : p->getInteriorRing(j-1);
                for (size_t k=0; k < pLinearRing->getNumPoints(); k++) {
                    x =  pLinearRing->getX(k);
                    y =  pLinearRing->getY(k);
                    pc->points.push_back(gda::Point(x,y));
                    main_map->set_bbox(x, y);
                }
            }
            main_map->records.push_back(pc);

        } else if (eType == wkbMultiPolygon) {
            gda::PolygonContents* pc = new gda::PolygonContents();
            OGRMultiPolygon* mpolygon = (OGRMultiPolygon*) geometry;

            int n_geom = mpolygon->getNumGeometries();

            // if there is more than one polygon, then we need to count
            // which part is processing accumulatively
            int numPoints = 0;
            OGRLinearRing* pLinearRing = 0;

            for (size_t i = 0; i < n_geom; i++ ) {
                OGRGeometry *ogrGeom = mpolygon->getGeometryRef(i);
                OGRPolygon *p = (OGRPolygon *) ogrGeom;
                if (i == 0) {
                    CopyEnvelope(p, pc);
                } else {
                    OGREnvelope pBox;
                    p->getEnvelope(&pBox);
                    if (pc->box[0] > pBox.MinX) pc->box[0] = pBox.MinX;
                    if (pc->box[1] > pBox.MinY) pc->box[1] = pBox.MinY;
                    if (pc->box[2] < pBox.MaxX) pc->box[2] = pBox.MaxX;
                    if (pc->box[3] < pBox.MaxY) pc->box[3] = pBox.MaxY;
                }
                // number of interior rings + 1 exterior ring
                int ni_rings = p->getNumInteriorRings() + 1;
                pc->num_parts += ni_rings;

                for (size_t j = 0; j < ni_rings; j++) {
                    pLinearRing = j == 0 ? p->getExteriorRing() : p->getInteriorRing(j - 1);
                    pc->holes.push_back(j>0);
                    pc->parts.push_back(numPoints);
                    numPoints += pLinearRing->getNumPoints();
                }
            }

            // resize points memory
            pc->num_points = numPoints;

            int pidx =0;
            double x, y;

            for (size_t i = 0; i < n_geom; i++ ) {
                OGRGeometry *ogrGeom = mpolygon->getGeometryRef(i);
                OGRPolygon *p = (OGRPolygon *) ogrGeom;

                // number of interior rings + 1 exterior ring
                int ni_rings = p->getNumInteriorRings() + 1;

                // read points
                for (size_t j=0; j < ni_rings; j++) {
                    pLinearRing = j==0 ? p->getExteriorRing() : p->getInteriorRing(j-1);
                    for (int k=0; k < pLinearRing->getNumPoints(); k++) {
                        x = pLinearRing->getX(k);
                        y = pLinearRing->getY(k);
                        main_map->set_bbox(x,y);
                        pc->points.push_back(gda::Point(x,y));
                    }
                }
            }
            main_map->records.push_back(pc);

        } else {
            std::string open_err_msg = "GeoDa does not support datasource with line data at this time.  Please choose"
                                       " a datasource with either point or polygon data.";
        }
    }

    return has_null_geometry;
     */
}
