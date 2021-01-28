#ifndef __POST_GEOMS_H__
#define __POST_GEOMS_H__

#ifdef _MSC_VER
typedef signed __int8 int8_t;
typedef unsigned __int8 uint8_t;
typedef signed __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

#if (_MSC_VER == 1500) || defined(__PG_GEODSA__)
//typedef unsigned __int32 size_t;
#define inline __inline
#define isnan _isnan
#endif



#define NULL 0
#else
#include <stdint.h>
#include <stdlib.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

/**
* LWTYPE numbers, used internally by PostGIS
*/
#define POINTTYPE                1
#define LINETYPE                 2
#define POLYGONTYPE              3
#define MULTIPOINTTYPE           4
#define MULTILINETYPE            5
#define MULTIPOLYGONTYPE         6
#define COLLECTIONTYPE           7
#define CIRCSTRINGTYPE           8
#define COMPOUNDTYPE             9
#define CURVEPOLYTYPE           10
#define MULTICURVETYPE          11
#define MULTISURFACETYPE        12
#define POLYHEDRALSURFACETYPE   13
#define TRIANGLETYPE            14
#define TINTYPE                 15
#define NUMTYPES                16

/**
* Macros for manipulating the 'flags' byte. A uint8_t used as follows:
* VVSRGBMZ
* Version bit, followed by
* Validty, Solid, ReadOnly, Geodetic, HasBBox, HasM and HasZ flags.
*/
#define LWFLAG_Z        0x01
#define LWFLAG_M        0x02
#define LWFLAG_BBOX     0x04
#define LWFLAG_GEODETIC 0x08
#define LWFLAG_READONLY 0x10
#define LWFLAG_SOLID    0x20

#define FLAGS_GET_Z(flags)         ((flags) & LWFLAG_Z)
#define FLAGS_GET_M(flags)        (((flags) & LWFLAG_M)>>1)
#define FLAGS_GET_BBOX(flags)     (((flags) & LWFLAG_BBOX)>>2)
#define FLAGS_GET_GEODETIC(flags) (((flags) & LWFLAG_GEODETIC)>>3)
#define FLAGS_GET_READONLY(flags) (((flags) & LWFLAG_READONLY)>>4)
#define FLAGS_GET_SOLID(flags)    (((flags) & LWFLAG_SOLID)>>5)

#define FLAGS_SET_Z(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_Z) : ((flags) & ~LWFLAG_Z))
#define FLAGS_SET_M(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_M) : ((flags) & ~LWFLAG_M))
#define FLAGS_SET_BBOX(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_BBOX) : ((flags) & ~LWFLAG_BBOX))
#define FLAGS_SET_GEODETIC(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_GEODETIC) : ((flags) & ~LWFLAG_GEODETIC))
#define FLAGS_SET_READONLY(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_READONLY) : ((flags) & ~LWFLAG_READONLY))
#define FLAGS_SET_SOLID(flags, value) ((flags) = (value) ? ((flags) | LWFLAG_SOLID) : ((flags) & ~LWFLAG_SOLID))

#define FLAGS_NDIMS(flags) (2 + FLAGS_GET_Z(flags) + FLAGS_GET_M(flags))
#define FLAGS_GET_ZM(flags) (FLAGS_GET_M(flags) + FLAGS_GET_Z(flags) * 2)
#define FLAGS_NDIMS_BOX(flags) (FLAGS_GET_GEODETIC(flags) ? 3 : FLAGS_NDIMS(flags))

/**
* Macro for reading the size from the GSERIALIZED size attribute.
* Cribbed from PgSQL, top 30 bits are size. Use VARSIZE() when working
* internally with PgSQL. See SET_VARSIZE_4B / VARSIZE_4B in
* PGSRC/src/include/postgres.h for details.
*/
#ifdef WORDS_BIGENDIAN
#define SIZE_GET(varsize) ((varsize) & 0x3FFFFFFF)
#define SIZE_SET(varsize, len) ((varsize) = ((len) & 0x3FFFFFFF))
#define IS_BIG_ENDIAN 1
#else
#define SIZE_GET(varsize) (((varsize) >> 2) & 0x3FFFFFFF)
#define SIZE_SET(varsize, len) ((varsize) = (((uint32_t)(len)) << 2))
#define IS_BIG_ENDIAN 0
#endif

/**
 * Parser check flags
 *
 *  @see lwgeom_from_wkb
 *  @see lwgeom_from_hexwkb
 *  @see lwgeom_parse_wkt
 */
#define LW_PARSER_CHECK_MINPOINTS  1
#define LW_PARSER_CHECK_ODD        2
#define LW_PARSER_CHECK_CLOSURE    4
#define LW_PARSER_CHECK_ZCLOSURE   8

#define LW_PARSER_CHECK_NONE   0
#define LW_PARSER_CHECK_ALL (LW_PARSER_CHECK_MINPOINTS | LW_PARSER_CHECK_ODD | LW_PARSER_CHECK_CLOSURE)

/**
* Return types for functions with status returns.
*/
#define LW_TRUE 1
#define LW_FALSE 0
#define LW_UNKNOWN 2
#define LW_FAILURE 0
#define LW_SUCCESS 1


/*
* this will change to NaN when I figure out how to
* get NaN in a platform-independent way
*/
#define NO_VALUE 0.0
#define NO_Z_VALUE NO_VALUE
#define NO_M_VALUE NO_VALUE

/**
* Flags applied in EWKB to indicate Z/M dimensions and
* presence/absence of SRID and bounding boxes
*/
#define WKBZOFFSET  0x80000000
#define WKBMOFFSET  0x40000000
#define WKBSRIDFLAG 0x20000000
#define WKBBBOXFLAG 0x10000000


/**
* Well-Known Binary (WKB) Output Variant Types
*/

#define WKB_DOUBLE_SIZE 8 /* Internal use only */
#define WKB_INT_SIZE 4 /* Internal use only */
#define WKB_BYTE_SIZE 1 /* Internal use only */

/**
* Well-Known Binary (WKB) Geometry Types
*/
#define WKB_POINT_TYPE 1
#define WKB_LINESTRING_TYPE 2
#define WKB_POLYGON_TYPE 3
#define WKB_MULTIPOINT_TYPE 4
#define WKB_MULTILINESTRING_TYPE 5
#define WKB_MULTIPOLYGON_TYPE 6
#define WKB_GEOMETRYCOLLECTION_TYPE 7
#define WKB_CIRCULARSTRING_TYPE 8
#define WKB_COMPOUNDCURVE_TYPE 9
#define WKB_CURVEPOLYGON_TYPE 10
#define WKB_MULTICURVE_TYPE 11
#define WKB_MULTISURFACE_TYPE 12
#define WKB_CURVE_TYPE 13 /* from ISO draft, not sure is real */
#define WKB_SURFACE_TYPE 14 /* from ISO draft, not sure is real */
#define WKB_POLYHEDRALSURFACE_TYPE 15
#define WKB_TIN_TYPE 16
#define WKB_TRIANGLE_TYPE 17


/* Memory management */
extern void *lwalloc(size_t size);
extern void *lwrealloc(void *mem, size_t size);
extern void lwfree(void *mem);

/**
* Used for passing the parse state between the parsing functions.
*/
typedef struct
{
    const uint8_t *wkb; /* Points to start of WKB */
    int32_t srid;    /* Current SRID we are handling */
    size_t wkb_size; /* Expected size of WKB */
    int8_t swap_bytes;  /* Do an endian flip? */
    int8_t check;       /* Simple validity checks on geometries */
    int8_t lwtype;      /* Current type we are handling */
    int8_t has_z;       /* Z? */
    int8_t has_m;       /* M? */
    int8_t has_srid;    /* SRID? */
    int8_t error;       /* An error was found (not enough bytes to read) */
    const uint8_t *pos; /* Current parse position */
} wkb_parse_state;

/**
* Check that we are not about to read off the end of the WKB
* array.
*/
static inline void wkb_parse_state_check(wkb_parse_state *s, size_t next)
{
    if( (s->pos + next) > (s->wkb + s->wkb_size) )
    {
        //lwerror("WKB structure does not match expected size!");
        s->error = LW_TRUE;
    }
}

/******************************************************************
* LWGEOM and GBOX both use LWFLAGS bit mask.
* Serializations (may) use different bit mask schemes.
*/
typedef uint16_t lwflags_t;

/******************************************************************
* GBOX structure.
* We include the flags (information about dimensionality),
* so we don't have to constantly pass them
* into functions that use the GBOX.
*/
typedef struct
{
    lwflags_t flags;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    double mmin;
    double mmax;
} GBOX;

/******************************************************************
* POINT2D, POINT3D, POINT3DM, POINT4D
*/
typedef struct
{
    double x, y;
}
        POINT2D;

typedef struct
{
    double x, y, z;
}
        POINT3DZ;

typedef struct
{
    double x, y, z;
}
        POINT3D;

typedef struct
{
    double x, y, m;
}
        POINT3DM;

typedef struct
{
    double x, y, z, m;
}
        POINT4D;

/******************************************************************
*  POINTARRAY
*  Point array abstracts a lot of the complexity of points and point lists.
*  It handles 2d/3d translation
*    (2d points converted to 3d will have z=0 or NaN)
*  DO NOT MIX 2D and 3D POINTS! EVERYTHING* is either one or the other
*/
typedef struct
{
    uint32_t npoints;   /* how many points we are currently storing */
    uint32_t maxpoints; /* how many points we have space for in serialized_pointlist */

    /* Use FLAGS_* macros to handle */
    lwflags_t flags;

    /* Array of POINT 2D, 3D or 4D, possibly misaligned. */
    uint8_t *serialized_pointlist;
}
        POINTARRAY;

/******************************************************************
* GSERIALIZED
*/

typedef struct
{
    uint32_t size; /* For PgSQL use only, use VAR* macros to manipulate. */
    uint8_t srid[3]; /* 24 bits of SRID */
    uint8_t gflags; /* HasZ, HasM, HasBBox, IsGeodetic */
    uint8_t data[1]; /* See gserialized.txt */
} GSERIALIZED;

/******************************************************************
* LWGEOM (any geometry type)
*
* Abstract type, note that 'type', 'bbox' and 'srid' are available in
* all geometry variants.
*/
typedef struct
{
    GBOX *bbox;
    void *data;
    int32_t srid;
    lwflags_t flags;
    uint8_t type;
    char pad[1]; /* Padding to 24 bytes (unused) */
}
        LWGEOM;

/* POINTYPE */
typedef struct
{
    GBOX *bbox;
    POINTARRAY *point;  /* hide 2d/3d (this will be an array of 1 point) */
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* POINTTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
}
LWPOINT; /* "light-weight point" */

 /* POLYGONTYPE */
typedef struct
{
    GBOX *bbox;
    POINTARRAY **rings; /* list of rings (list of points) */
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* POLYGONTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
    uint32_t nrings;   /* how many rings we are currently storing */
    uint32_t maxrings; /* how many rings we have space for in **rings */
}
        LWPOLY; /* "light-weight polygon" */

/* MULTIPOINTTYPE */
typedef struct
{
    GBOX *bbox;
    LWPOINT **geoms;
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* MULTYPOINTTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
    uint32_t ngeoms;   /* how many geometries we are currently storing */
    uint32_t maxgeoms; /* how many geometries we have space for in **geoms */
}
        LWMPOINT;

/* LINETYPE */
typedef struct
{
    GBOX *bbox;
    POINTARRAY *points; /* array of POINT3D */
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* LINETYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
}
        LWLINE; /* "light-weight line" */

/* MULTIPOLYGONTYPE */
typedef struct
{
    GBOX *bbox;
    LWPOLY **geoms;
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* MULTIPOLYGONTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
    uint32_t ngeoms;   /* how many geometries we are currently storing */
    uint32_t maxgeoms; /* how many geometries we have space for in **geoms */
}
        LWMPOLY;

/* COLLECTIONTYPE */
typedef struct
{
    GBOX *bbox;
    LWGEOM **geoms;
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* COLLECTIONTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
    uint32_t ngeoms;   /* how many geometries we are currently storing */
    uint32_t maxgeoms; /* how many geometries we have space for in **geoms */
}
        LWCOLLECTION;

/* CURVEPOLYTYPE */
typedef struct
{
    GBOX *bbox;
    LWGEOM **rings;
    int32_t srid;
    lwflags_t flags;
    uint8_t type; /* CURVEPOLYTYPE */
    char pad[1]; /* Padding to 24 bytes (unused) */
    uint32_t nrings;    /* how many rings we are currently storing */
    uint32_t maxrings;  /* how many rings we have space for in **rings */
}
        LWCURVEPOLY; /* "light-weight polygon" */

/**
* Create a new gbox with the dimensionality indicated by the flags. Caller
* is responsible for freeing.
*/
extern GBOX* gbox_new(lwflags_t flags);

/*
 * Size of point represeneted in the POINTARRAY
 * 16 for 2d, 24 for 3d, 32 for 4d
 */
static inline size_t
ptarray_point_size(const POINTARRAY *pa)
{
    return sizeof(double) * FLAGS_NDIMS(pa->flags);
}

/****************************************************************
 * MEMORY MANAGEMENT
 ****************************************************************/

/*
* The *_free family of functions frees *all* memory associated
* with the pointer. When the recursion gets to the level of the
* POINTARRAY, the POINTARRAY is only freed if it is not flagged
* as "read only". LWGEOMs constructed on top of GSERIALIZED
* from PgSQL use read only point arrays.
*/

extern void ptarray_free(POINTARRAY *pa);
extern void lwpoint_free(LWPOINT *pt);
extern void lwpoly_free(LWPOLY *poly);
extern void lwmpoint_free(LWMPOINT *mpt);
extern void lwmpoly_free(LWMPOLY *mpoly);
extern void lwcollection_free(LWCOLLECTION *col);
extern void lwgeom_free(LWGEOM *geom);

/**
* Construct a new #POINTARRAY, <em>copying</em> in the data from ptlist
*/
extern POINTARRAY* ptarray_construct_copy_data(char hasz, char hasm, uint32_t npoints, const uint8_t *ptlist);

/**
* Construct an empty pointarray, allocating storage and setting
* the npoints, but not filling in any information. Should be used in conjunction
* with ptarray_set_point4d to fill in the information in the array.
*/
extern POINTARRAY* ptarray_construct(char hasz, char hasm, uint32_t npoints);

/**
* Create a new #POINTARRAY with no points. Allocate enough storage
* to hold maxpoints vertices before having to reallocate the storage
* area.
*/
extern POINTARRAY* ptarray_construct_empty(char hasz, char hasm, uint32_t maxpoints);

/*
* Empty geometry constructors.
*/
//extern LWPOINT* lwpoint_construct_empty(int32_t srid, char hasz, char hasm);

/*
* Geometry constructors. These constructors to not copy the point arrays
* passed to them, they just take references, so do not free them out
* from underneath the geometries.
*/
extern LWPOINT* lwpoint_construct(int32_t srid, GBOX *bbox, POINTARRAY *point);

/**
* Return the type name string associated with a type number
* (e.g. Point, LineString, Polygon)
*/
extern const char *lwtype_name(uint8_t type);

/*
 * Get a pointer to Nth point of a POINTARRAY
 * You'll need to cast it to appropriate dimensioned point.
 * Note that if you cast to a higher dimensional point you'll
 * possibly corrupt the POINTARRAY.
 *
 * Casting to returned pointer to POINT2D* should be safe,
 * as gserialized format always keeps the POINTARRAY pointer
 * aligned to double boundary.
 *
 * WARNING: Don't cast this to a POINT!
 * it would not be reliable due to memory alignment constraints
 */
static inline uint8_t *
getPoint_internal(const POINTARRAY *pa, uint32_t n)
{
    size_t size;
    uint8_t *ptr;

#if PARANOIA_LEVEL > 0
    assert(pa);
	assert(n <= pa->npoints);
	assert(n <= pa->maxpoints);
#endif

    size = ptarray_point_size(pa);
    ptr = pa->serialized_pointlist + size * n;

    return ptr;
}
/**
 * Returns a POINT2D pointer into the POINTARRAY serialized_ptlist,
 * suitable for reading from. This is very high performance
 * and declared const because you aren't allowed to muck with the
 * values, only read them.
 */
static inline const POINT2D *
getPoint2d_cp(const POINTARRAY *pa, uint32_t n)
{
    return (const POINT2D *)getPoint_internal(pa, n);
}

/**
* Add a ring, allocating extra space if necessary. The polygon takes
* ownership of the passed point array.
*/
extern int lwpoly_add_ring(LWPOLY *poly, POINTARRAY *pa);

/**
 * @param wkb_size length of WKB byte buffer
 * @param wkb WKB byte buffer
 * @param check parser check flags, see LW_PARSER_CHECK_* macros
 */
extern LWGEOM* lwgeom_from_wkb(const uint8_t *wkb, const size_t wkb_size, const char check);

/*
* Empty geometry constructors.
*/
//extern LWGEOM* lwgeom_construct_empty(uint8_t type, int32_t srid, char hasz, char hasm);
extern LWPOINT* lwpoint_construct_empty(int32_t srid, char hasz, char hasm);
extern LWPOLY* lwpoly_construct_empty(int32_t srid, char hasz, char hasm);
extern LWCURVEPOLY* lwcurvepoly_construct_empty(int32_t srid, char hasz, char hasm);
//extern LWMPOINT* lwmpoint_construct_empty(int32_t srid, char hasz, char hasm);
//extern LWMPOLY* lwmpoly_construct_empty(int32_t srid, char hasz, char hasm);
extern LWCOLLECTION* lwcollection_construct_empty(uint8_t type, int32_t srid, char hasz, char hasm);

/* Casts LWGEOM->LW* (return NULL if cast is illegal) */
static inline LWPOINT *
lwgeom_as_lwpoint(const LWGEOM *lwgeom)
{
    if (!lwgeom)
        return NULL;
    if (lwgeom->type == POINTTYPE)
        return (LWPOINT *)lwgeom;
    else
        return NULL;
}

LWMPOLY *lwgeom_as_lwmpoly(const LWGEOM *lwgeom);
LWMPOINT *lwgeom_as_lwmpoint(const LWGEOM *lwgeom);
LWPOLY *lwgeom_as_lwpoly(const LWGEOM *lwgeom);

/******************************************************************/
/* Functions that work on type numbers */

/**
* Determine whether a type number is a collection or not
*/
extern int lwtype_is_collection(uint8_t type);

extern LWCOLLECTION* lwcollection_add_lwgeom(LWCOLLECTION *col, const LWGEOM *geom);

/** Check if subtype is allowed in collectiontype */
int lwcollection_allows_subtype(int collectiontype, int subtype);

/** Ensure the collection can hold at least up to ngeoms geometries */
void lwcollection_reserve(LWCOLLECTION *col, uint32_t ngeoms);

/**
* Return a valid SRID from an arbitrary integer
* Raises a notice if what comes out is different from
* what went in.
* Raises an error if SRID value is out of bounds.
*/
extern int32_t clamp_srid(int32_t srid);


/**
* Deep clone an LWGEOM, everything is copied
*/
extern LWGEOM *lwgeom_clone_deep(const LWGEOM *lwgeom);
extern POINTARRAY *ptarray_clone_deep(const POINTARRAY *ptarray);
extern int ptarray_is_closed_2d(const POINTARRAY *pa);

LWLINE *lwline_clone_deep(const LWLINE *lwgeom);
LWPOLY *lwpoly_clone_deep(const LWPOLY *lwgeom);
LWCOLLECTION *lwcollection_clone_deep(const LWCOLLECTION *lwgeom);

/**
* Return a copy of the #GBOX, based on dimensionality of flags.
*/
 GBOX* gbox_copy(const GBOX *gbox);

/*
 * copies a point from the point array into the parameter point
 * will set point's z=0 (or NaN) if pa is 2d
 * will set point's m=0 (or NaN) if pa is 3d or 2d
 * NOTE: point is a real POINT3D *not* a pointer
 */
POINT4D getPoint4d(const POINTARRAY *pa, uint32_t n);
/*
 * copies a point from the point array into the parameter point
 * will set point's z=0 (or NaN) if pa is 2d
 * will set point's m=0 (or NaN) if pa is 3d or 2d
 * NOTE: this will modify the point4d pointed to by 'point'.
 */
 int getPoint4d_p(const POINTARRAY *pa, uint32_t n, POINT4D *point);

/**
 * Return true or false depending on whether a geometry is an "empty"
 * geometry (no vertices members)
 */
static inline int
lwpoint_is_empty(const LWPOINT *point)
{
    return !point->point || point->point->npoints < 1;
}

static inline int
lwpoly_is_empty(const LWPOLY *poly)
{
    return poly->nrings < 1 || !poly->rings || !poly->rings[0] || poly->rings[0]->npoints < 1;
}

static inline int lwgeom_is_empty(const LWGEOM *geom);


static inline int
lwcollection_is_empty(const LWCOLLECTION *col)
{
    uint32_t i;
    if (col->ngeoms == 0 || !col->geoms)
        return LW_TRUE;
    for (i = 0; i < col->ngeoms; i++)
    {
        if (!lwgeom_is_empty(col->geoms[i]))
            return LW_FALSE;
    }
    return LW_TRUE;
}



static inline int
lwgeom_is_empty(const LWGEOM *geom)
{
    switch (geom->type)
    {
        case POINTTYPE:
            return lwpoint_is_empty((LWPOINT *)geom);
            break;
        case LINETYPE:
            return LW_FALSE;
            break;
        case CIRCSTRINGTYPE:
            return LW_FALSE;//lwcircstring_is_empty((LWCIRCSTRING *)geom);
            break;
        case POLYGONTYPE:
            return lwpoly_is_empty((LWPOLY *)geom);
            break;
        case TRIANGLETYPE:
            return LW_FALSE;//lwtriangle_is_empty((LWTRIANGLE *)geom);
            break;
        case MULTIPOINTTYPE:
        case MULTILINETYPE:
        case MULTIPOLYGONTYPE:
        case COMPOUNDTYPE:
        case CURVEPOLYTYPE:
        case MULTICURVETYPE:
        case MULTISURFACETYPE:
        case POLYHEDRALSURFACETYPE:
        case TINTYPE:
        case COLLECTIONTYPE:
            return lwcollection_is_empty((LWCOLLECTION *)geom);
            break;
        default:
            return LW_FALSE;
            break;
    }
}


#ifdef __cplusplus
}
#endif

#endif
