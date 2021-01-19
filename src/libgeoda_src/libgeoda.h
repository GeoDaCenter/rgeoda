#ifndef LIBGEODA_LIBRARY_H
#define LIBGEODA_LIBRARY_H

#ifdef __MSVC__
#include <Windows.h>
#include <Psapi.h>
#endif

#include <list>
#include <vector>
#include <algorithm>
#include <map>

#include "./pg/geoms.h"
#include "gda_interface.h"

// forward declaration



class GeoDaColumn {
public:
    enum FieldType { integer_type, string_type, real_type };
    std::string name;
    FieldType field_type;
    int field_length;
    int field_decimals;
    std::vector<bool> undefs;

    GeoDaColumn(const std::string& nm, FieldType ft, int flen, int fdec)
    : name(nm), field_type(ft), field_length(flen), field_decimals(fdec) {}
    virtual ~GeoDaColumn() {}
};

class GeoDaIntColumn : public GeoDaColumn {
public:
    std::vector<long long> data;

    virtual std::vector<long long>& GetData() { return data;}
    virtual void SetData(const std::vector<long long>& vals) { data = vals; }

    GeoDaIntColumn(const std::string& nm, const std::vector<long long>& vals)
    : GeoDaColumn(nm, integer_type, 20, 0), data(vals) {}
};


class GeoDaStringColumn : public GeoDaColumn {
public:
    std::vector<std::string> data;

    virtual std::vector<std::string>& GetData() { return data;}
    virtual void SetData(const std::vector<std::string>& vals) { data = vals; }

    GeoDaStringColumn(const std::string& nm, const std::vector<std::string>& vals)
            : GeoDaColumn(nm, string_type, 254, 0), data(vals) {}
};

class GeoDaRealColumn : public GeoDaColumn {
public:
    std::vector<double> data;

    virtual std::vector<double>& GetData() { return data;}
    virtual void SetData(const std::vector<double>& vals) { data = vals; }

    GeoDaRealColumn(const std::string& nm, const std::vector<double>& vals)
            : GeoDaColumn(nm, real_type, 35, 15), data(vals) {}
};

class GeoDaTable {
public:
    GeoDaTable(){};
    virtual ~GeoDaTable(){};

    void AddIntColumn(const std::string& nm,
            const std::vector<long long>& vals,
            const std::vector<bool>& undefs=std::vector<bool>()) {
        GeoDaColumn* col = new GeoDaIntColumn(nm, vals);
        col->undefs = undefs;
        columns.push_back(col);
    }
    void AddStringColumn(const std::string& nm,
            const std::vector<std::string>& vals,
            const std::vector<bool>& undefs=std::vector<bool>()) {
        GeoDaColumn* col = new GeoDaStringColumn(nm, vals);
        col->undefs = undefs;
        columns.push_back(col);
    }
    void AddRealColumn(const std::string& nm,
            const std::vector<double>& vals,
            const std::vector<bool>& undefs=std::vector<bool>()) {
        GeoDaColumn* col = new GeoDaRealColumn(nm, vals);
        col->undefs = undefs;
        columns.push_back(col);
    }

    GeoDaColumn* GetColumn(const std::string& col_name) {
        for (size_t i=0; i<columns.size(); ++i) {
            if (col_name.compare(columns[i]->name) == 0) {
                return columns[i];
            }
        }
        return 0;
    }

    GeoDaColumn* GetColumn(int idx) {return columns[idx];}

    int GetNumCols() { return columns.size(); }

protected:
    std::vector<GeoDaColumn*> columns;
};

class GeoDa : public AbstractGeoDa {
public:
    enum MapType { point_type, polygon_type, line_type, unknown_type };

    // this constructor is for Python
    GeoDa(const std::string& layer_name,
          const std::string& map_type,
          const std::vector<unsigned char> &wkbs,
          const std::vector<int>& wkb_bytes_len);

	// this constructor is for R
    GeoDa(const std::string& layer_name,
          const std::string& map_type,
          int num_features,
          unsigned char* wkbs,
          const std::vector<int>& wkb_bytes_len);

    // this constructor is for ESRI Shapefile
    GeoDa(const char* pDsPath, const char* layer_name=NULL);

    virtual ~GeoDa();

    // interfaces from AbstractGeoDa
    virtual int GetNumObs() const;
    virtual const std::vector<gda::PointContents*>& GetCentroids();
    virtual int GetMapType();
    virtual gda::MainMap& GetMainMap();

    // Layer functions
    int GetNumCols() const;
    std::vector<std::string> GetFieldTypes();
    std::vector<std::string> GetFieldNames();

    // Data functions
    std::vector<double> GetNumericCol(std::string col_name);
    std::vector<long long> GetIntegerCol(std::string col_name);
    std::vector<std::string> GetStringCol(std::string col_name);
    std::vector<bool> GetNullValues(std::string col_name);


protected:
    void AddPoint(LWPOINT* lw_pt);
    void AddMultiPoint( LWMPOINT* lw_mpt);
    void AddPolygon( LWPOLY* lw_poly);
    void AddMultiPolygon( LWMPOLY* lw_mpoly);
    void AddNullGeometry();

    void ReadShapefile(const char* fpath);
    void ReadDbffile(const char* fpath);

    // reserved for working with GDAL
    void ReadAllFeatures();

	void Init(const std::string& layer_name,
            const std::string& map_type,
            int num_features,
            unsigned char* wkbs,
            const std::vector<int>& wkb_bytes_len);

    static const std::string DT_STRING;
    static const std::string DT_INTEGER;
    static const std::string DT_NUMERIC;

    MapType mapType;

    int numLayers;
    int numObs;
    int numCols;

    GeoDaTable* table;
    std::vector<std::string> fieldNames;
    std::vector<std::string> fieldTypes;
    std::map<std::string, unsigned int> fieldNameIdx;

    std::vector<gda::PointContents*> centroids;

    gda::MainMap* main_map;
};

int test();

GeoDaColumn* ToGeoDaColumn(GeoDaStringColumn* col);
GeoDaColumn* ToGeoDaColumn(GeoDaIntColumn* col);
GeoDaColumn* ToGeoDaColumn(GeoDaRealColumn* col);


GeoDa* CreateGeoDaFromGPD(const std::string& layer_name,
                         const std::string& map_type,
                         const std::vector<unsigned char> &wkbs,
                         const std::vector<int>& wkb_bytes_len);

GeoDa* CreateGeoDaFromSHP(const char* pDsPath, const char* layer_name=NULL);

#endif
