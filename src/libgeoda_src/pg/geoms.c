#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "config.h"
#include "utils.h"
#include "geoms.h"


/**
* Internal function declarations.
*/
LWGEOM* lwgeom_from_wkb_state(wkb_parse_state *s);

static char *lwgeomTypeName[] =
        {
                "Unknown",
                "Point",
                "LineString",
                "Polygon",
                "MultiPoint",
                "MultiLineString",
                "MultiPolygon",
                "GeometryCollection",
                "CircularString",
                "CompoundCurve",
                "CurvePolygon",
                "MultiCurve",
                "MultiSurface",
                "PolyhedralSurface",
                "Triangle",
                "Tin"
        };

/* Default allocators */
static void * default_allocator(size_t size);
static void default_freeor(void *mem);
static void * default_reallocator(void *mem, size_t size);
lwallocator lwalloc_var = default_allocator;
lwreallocator lwrealloc_var = default_reallocator;
lwfreeor lwfree_var = default_freeor;

/*
 * Default allocators
 *
 * We include some default allocators that use malloc/free/realloc
 * along with stdout/stderr since this is the most common use case
 *
 */

static void
default_freeor(void *mem)
{
    free(mem);
}

static void *
default_allocator(size_t size)
{
    void *mem = malloc(size);
    return mem;
}

static void *
default_reallocator(void *mem, size_t size)
{
    void *ret = realloc(mem, size);
    return ret;
}

void *
lwalloc(size_t size)
{
    void *mem = lwalloc_var(size);
    LWDEBUGF(4, "lwalloc: %d@%p", size, mem);
    return mem;
}

void *
lwrealloc(void *mem, size_t size)
{
    LWDEBUGF(5, "lwrealloc: %d@%p", size, mem);
    return lwrealloc_var(mem, size);
}

void
lwfree(void *mem)
{
    lwfree_var(mem);
}



int32_t
clamp_srid(int32_t srid)
{
    int newsrid = srid;

    if ( newsrid <= 0 ) {
        if ( newsrid != SRID_UNKNOWN ) {
            newsrid = SRID_UNKNOWN;
            lwnotice("SRID value %d converted to the officially unknown SRID value %d", srid, newsrid);
        }
    } else if ( srid > SRID_MAXIMUM ) {
        newsrid = SRID_USER_MAXIMUM + 1 +
                  /* -1 is to reduce likelyhood of clashes */
                  /* NOTE: must match implementation in postgis_restore.pl */
                  ( srid % ( SRID_MAXIMUM - SRID_USER_MAXIMUM - 1 ) );
        lwnotice("SRID value %d > SRID_MAXIMUM converted to %d", srid, newsrid);
    }

    return newsrid;
}


lwflags_t lwflags(int hasz, int hasm, int geodetic)
{
    lwflags_t flags = 0;
    if (hasz)
        FLAGS_SET_Z(flags, 1);
    if (hasm)
        FLAGS_SET_M(flags, 1);
    if (geodetic)
        FLAGS_SET_GEODETIC(flags, 1);
    return flags;
}

const char*
lwtype_name(uint8_t type)
{
    if ( type > 15 )
    {
        /* assert(0); */
        return "Invalid type";
    }
    return lwgeomTypeName[(int ) type];
}

/** Return TRUE if the geometry may contain sub-geometries, i.e. it is a MULTI* or COMPOUNDCURVE */
int
lwtype_is_collection(uint8_t type)
{

    switch (type)
    {
        case MULTIPOINTTYPE:
        case MULTILINETYPE:
        case MULTIPOLYGONTYPE:
        case COLLECTIONTYPE:
        case CURVEPOLYTYPE:
        case COMPOUNDTYPE:
        case MULTICURVETYPE:
        case MULTISURFACETYPE:
        case POLYHEDRALSURFACETYPE:
        case TINTYPE:
            return LW_TRUE;
            break;

        default:
            return LW_FALSE;
    }
}

void
ptarray_free(POINTARRAY *pa)
{
    if (pa)
    {
        if (pa->serialized_pointlist && (!FLAGS_GET_READONLY(pa->flags)))
            lwfree(pa->serialized_pointlist);
        lwfree(pa);
    }
}

void lwpoint_free(LWPOINT *pt)
{
    if ( ! pt ) return;

    if ( pt->bbox )
        lwfree(pt->bbox);
    if ( pt->point )
        ptarray_free(pt->point);
    lwfree(pt);
}

LWPOLY *
lwpoly_construct_empty(int32_t srid, char hasz, char hasm)
{
    LWPOLY *result = lwalloc(sizeof(LWPOLY));
    result->type = POLYGONTYPE;
    result->flags = lwflags(hasz,hasm,0);
    result->srid = srid;
    result->nrings = 0;
    result->maxrings = 1; /* Allocate room for ring, just in case. */
    result->rings = lwalloc(result->maxrings * sizeof(POINTARRAY*));
    result->bbox = NULL;
    return result;
}

int lwcollection_allows_subtype(int collectiontype, int subtype)
{
    if ( collectiontype == COLLECTIONTYPE )
        return LW_TRUE;
    if ( collectiontype == MULTIPOINTTYPE &&
         subtype == POINTTYPE )
        return LW_TRUE;
    if ( collectiontype == MULTILINETYPE &&
         subtype == LINETYPE )
        return LW_TRUE;
    if ( collectiontype == MULTIPOLYGONTYPE &&
         subtype == POLYGONTYPE )
        return LW_TRUE;
    if ( collectiontype == COMPOUNDTYPE &&
         (subtype == LINETYPE || subtype == CIRCSTRINGTYPE) )
        return LW_TRUE;
    if ( collectiontype == CURVEPOLYTYPE &&
         (subtype == CIRCSTRINGTYPE || subtype == LINETYPE || subtype == COMPOUNDTYPE) )
        return LW_TRUE;
    if ( collectiontype == MULTICURVETYPE &&
         (subtype == CIRCSTRINGTYPE || subtype == LINETYPE || subtype == COMPOUNDTYPE) )
        return LW_TRUE;
    if ( collectiontype == MULTISURFACETYPE &&
         (subtype == POLYGONTYPE || subtype == CURVEPOLYTYPE) )
        return LW_TRUE;
    if ( collectiontype == POLYHEDRALSURFACETYPE &&
         subtype == POLYGONTYPE )
        return LW_TRUE;
    if ( collectiontype == TINTYPE &&
         subtype == TRIANGLETYPE )
        return LW_TRUE;

    /* Must be a bad combination! */
    return LW_FALSE;
}

/**
 * Ensure the collection can hold up at least ngeoms
 */
void lwcollection_reserve(LWCOLLECTION *col, uint32_t ngeoms)
{
    if ( ngeoms <= col->maxgeoms ) return;

    /* Allocate more space if we need it */
    do { col->maxgeoms *= 2; } while ( col->maxgeoms < ngeoms );
    col->geoms = lwrealloc(col->geoms, sizeof(LWGEOM*) * col->maxgeoms);
}

void
lwpoly_free(LWPOLY* poly)
{
    uint32_t t;

    if (!poly) return;

    if (poly->bbox) lwfree(poly->bbox);

    if ( poly->rings )
    {
        for (t = 0; t < poly->nrings; t++)
            if (poly->rings[t]) ptarray_free(poly->rings[t]);
        lwfree(poly->rings);
    }

    lwfree(poly);
}

void lwmpoint_free(LWMPOINT *mpt)
{
    uint32_t i;

    if ( ! mpt ) return;

    if ( mpt->bbox )
        lwfree(mpt->bbox);

    for ( i = 0; i < mpt->ngeoms; i++ )
        if ( mpt->geoms && mpt->geoms[i] )
            lwpoint_free(mpt->geoms[i]);

    if ( mpt->geoms )
        lwfree(mpt->geoms);

    lwfree(mpt);
}

void lwmpoly_free(LWMPOLY *mpoly)
{
    uint32_t i;
    if ( ! mpoly ) return;
    if ( mpoly->bbox )
        lwfree(mpoly->bbox);

    for ( i = 0; i < mpoly->ngeoms; i++ )
        if ( mpoly->geoms && mpoly->geoms[i] )
            lwpoly_free(mpoly->geoms[i]);

    if ( mpoly->geoms )
        lwfree(mpoly->geoms);

    lwfree(mpoly);
}

void lwcollection_free(LWCOLLECTION *col)
{
    uint32_t i;
    if ( ! col ) return;

    if ( col->bbox )
    {
        lwfree(col->bbox);
    }
    for ( i = 0; i < col->ngeoms; i++ )
    {
        LWDEBUGF(4,"freeing geom[%d]", i);
        if ( col->geoms && col->geoms[i] )
            lwgeom_free(col->geoms[i]);
    }
    if ( col->geoms )
    {
        lwfree(col->geoms);
    }
    lwfree(col);
}

/**
* WKB inputs *must* have a declared size, to prevent malformed WKB from reading
* off the end of the memory segment (this stops a malevolent user from declaring
* a one-ring polygon to have 10 rings, causing the WKB reader to walk off the
* end of the memory).
*
* Check is a bitmask of: LW_PARSER_CHECK_MINPOINTS, LW_PARSER_CHECK_ODD,
* LW_PARSER_CHECK_CLOSURE, LW_PARSER_CHECK_NONE, LW_PARSER_CHECK_ALL
*/
LWGEOM* lwgeom_from_wkb(const uint8_t *wkb, const size_t wkb_size, const char check)
{
    wkb_parse_state s;

    /* Initialize the state appropriately */
    s.wkb = wkb;
    s.wkb_size = wkb_size;
    s.swap_bytes = LW_FALSE;
    s.check = check;
    s.lwtype = 0;
    s.srid = SRID_UNKNOWN;
    s.has_z = LW_FALSE;
    s.has_m = LW_FALSE;
    s.has_srid = LW_FALSE;
    s.error = LW_FALSE;
    s.pos = wkb;

    if (!wkb || !wkb_size)
        return NULL;

    return lwgeom_from_wkb_state(&s);
}

GBOX* gbox_new(lwflags_t flags)
{
    GBOX *g = (GBOX*)lwalloc(sizeof(GBOX));
    //gbox_init(g);
    memset(g, 0, sizeof(GBOX));
    g->flags = flags;
    return g;
}

LWMPOINT *
lwgeom_as_lwmpoint(const LWGEOM *lwgeom)
{
    if ( lwgeom == NULL ) return NULL;
    if ( lwgeom->type == MULTIPOINTTYPE )
        return (LWMPOINT *)lwgeom;
    else return NULL;
}


LWPOLY *
lwgeom_as_lwpoly(const LWGEOM *lwgeom)
{
    if ( lwgeom == NULL ) return NULL;
    if ( lwgeom->type == POLYGONTYPE )
        return (LWPOLY *)lwgeom;
    else return NULL;
}

LWMPOLY *
lwgeom_as_lwmpoly(const LWGEOM *lwgeom)
{
    if ( lwgeom == NULL ) return NULL;
    if ( lwgeom->type == MULTIPOLYGONTYPE )
        return (LWMPOLY *)lwgeom;
    else return NULL;
}



void lwgeom_free(LWGEOM *lwgeom)
{

    /* There's nothing here to free... */
    if( ! lwgeom ) return;

    LWDEBUGF(5,"freeing a %s",lwtype_name(lwgeom->type));

    switch (lwgeom->type)
    {
        case POINTTYPE:
            lwpoint_free((LWPOINT *)lwgeom);
            break;

        case POLYGONTYPE:
            lwpoly_free((LWPOLY *)lwgeom);
            break;

        case MULTIPOINTTYPE:
            lwmpoint_free((LWMPOINT *)lwgeom);
            break;

        case MULTIPOLYGONTYPE:
            lwmpoly_free((LWMPOLY *)lwgeom);
            break;

        case CURVEPOLYTYPE:
        case COMPOUNDTYPE:
        case MULTICURVETYPE:
        case MULTISURFACETYPE:
        case COLLECTIONTYPE:
            lwcollection_free((LWCOLLECTION *)lwgeom);
            break;
        case LINETYPE:
        case CIRCSTRINGTYPE:
        case TRIANGLETYPE:
        case MULTILINETYPE:
        case POLYHEDRALSURFACETYPE:
        case TINTYPE:
            lwerror("lwgeom_free called with unsupported type (%d) %s", lwgeom->type, lwtype_name(lwgeom->type));
        default:
            lwerror("lwgeom_free called with unknown type (%d) %s", lwgeom->type, lwtype_name(lwgeom->type));
    }
    return;
}

POINTARRAY*
ptarray_construct(char hasz, char hasm, uint32_t npoints)
{
    POINTARRAY *pa = ptarray_construct_empty(hasz, hasm, npoints);
    pa->npoints = npoints;
    return pa;
}

POINTARRAY*
ptarray_construct_empty(char hasz, char hasm, uint32_t maxpoints)
{
    POINTARRAY *pa = lwalloc(sizeof(POINTARRAY));
    pa->serialized_pointlist = NULL;

    /* Set our dimensionality info on the bitmap */
    pa->flags = lwflags(hasz, hasm, 0);

    /* We will be allocating a bit of room */
    pa->npoints = 0;
    pa->maxpoints = maxpoints;

    /* Allocate the coordinate array */
    if ( maxpoints > 0 )
        pa->serialized_pointlist = lwalloc(maxpoints * ptarray_point_size(pa));
    else
        pa->serialized_pointlist = NULL;

    return pa;
}

POINTARRAY*
ptarray_construct_copy_data(char hasz, char hasm, uint32_t npoints, const uint8_t *ptlist)
{
    POINTARRAY *pa = lwalloc(sizeof(POINTARRAY));

    pa->flags = lwflags(hasz, hasm, 0);
    pa->npoints = npoints;
    pa->maxpoints = npoints;

    if ( npoints > 0 )
    {
        pa->serialized_pointlist = lwalloc(ptarray_point_size(pa) * npoints);
        memcpy(pa->serialized_pointlist, ptlist, ptarray_point_size(pa) * npoints);
    }
    else
    {
        pa->serialized_pointlist = NULL;
    }

    return pa;
}

/*
 * Construct a new point.  point will not be copied
 * use SRID=SRID_UNKNOWN for unknown SRID (will have 8bit type's S = 0)
 */
LWPOINT *
lwpoint_construct(int32_t srid, GBOX *bbox, POINTARRAY *point)
{
    LWPOINT *result;
    lwflags_t flags = 0;

    if (point == NULL)
        return NULL; /* error */

    result = lwalloc(sizeof(LWPOINT));
    result->type = POINTTYPE;
    FLAGS_SET_Z(flags, FLAGS_GET_Z(point->flags));
    FLAGS_SET_M(flags, FLAGS_GET_M(point->flags));
    FLAGS_SET_BBOX(flags, bbox?1:0);
    result->flags = flags;
    result->srid = srid;
    result->point = point;
    result->bbox = bbox;

    return result;
}

LWPOINT *
lwpoint_construct_empty(int32_t srid, char hasz, char hasm)
{
    LWPOINT *result = lwalloc(sizeof(LWPOINT));
    result->type = POINTTYPE;
    result->flags = lwflags(hasz, hasm, 0);
    result->srid = srid;
    result->point = ptarray_construct(hasz, hasm, 0);
    result->bbox = NULL;
    return result;
}

LWCURVEPOLY *
lwcurvepoly_construct_empty(int32_t srid, char hasz, char hasm)
{
    LWCURVEPOLY *ret;

    ret = lwalloc(sizeof(LWCURVEPOLY));
    ret->type = CURVEPOLYTYPE;
    ret->flags = lwflags(hasz, hasm, 0);
    ret->srid = srid;
    ret->nrings = 0;
    ret->maxrings = 1; /* Allocate room for sub-members, just in case. */
    ret->rings = lwalloc(ret->maxrings * sizeof(LWGEOM*));
    ret->bbox = NULL;

    return ret;
}

int lwcurvepoly_add_ring(LWCURVEPOLY *poly, LWGEOM *ring)
{
    uint32_t i;

    /* Can't do anything with NULLs */
    if( ! poly || ! ring )
    {
        LWDEBUG(4,"NULL inputs!!! quitting");
        return LW_FAILURE;
    }

    /* Check that we're not working with garbage */
    if ( poly->rings == NULL && (poly->nrings || poly->maxrings) )
    {
        LWDEBUG(4,"mismatched nrings/maxrings");
        lwerror("Curvepolygon is in inconsistent state. Null memory but non-zero collection counts.");
        return LW_FAILURE;
    }

    /* Check that we're adding an allowed ring type */
    if ( ! ( ring->type == LINETYPE || ring->type == CIRCSTRINGTYPE || ring->type == COMPOUNDTYPE ) )
    {
        LWDEBUGF(4,"got incorrect ring type: %s",lwtype_name(ring->type));
        return LW_FAILURE;
    }


    /* In case this is a truly empty, make some initial space  */
    if ( poly->rings == NULL )
    {
        poly->maxrings = 2;
        poly->nrings = 0;
        poly->rings = lwalloc(poly->maxrings * sizeof(LWGEOM*));
    }

    /* Allocate more space if we need it */
    if ( poly->nrings == poly->maxrings )
    {
        poly->maxrings *= 2;
        poly->rings = lwrealloc(poly->rings, sizeof(LWGEOM*) * poly->maxrings);
    }

    /* Make sure we don't already have a reference to this geom */
    for ( i = 0; i < poly->nrings; i++ )
    {
        if ( poly->rings[i] == ring )
        {
            LWDEBUGF(4, "Found duplicate geometry in collection %p == %p", poly->rings[i], ring);
            return LW_SUCCESS;
        }
    }

    /* Add the ring and increment the ring count */
    poly->rings[poly->nrings] = (LWGEOM*)ring;
    poly->nrings++;
    return LW_SUCCESS;
}

LWCOLLECTION *
lwcollection_construct_empty(uint8_t type, int32_t srid, char hasz, char hasm)
{

    LWCOLLECTION* ret;
    if( ! lwtype_is_collection(type) )
    {
        lwerror("Non-collection type specified in collection constructor!");
        return NULL;
    }

    ret = lwalloc(sizeof(LWCOLLECTION));
    ret->type = type;
    ret->flags = lwflags(hasz,hasm,0);
    ret->srid = srid;
    ret->ngeoms = 0;
    ret->maxgeoms = 1; /* Allocate room for sub-members, just in case. */
    ret->geoms = lwalloc(ret->maxgeoms * sizeof(LWGEOM*));
    ret->bbox = NULL;

    return ret;
}

/**
* Byte
* Read a byte and advance the parse state forward.
*/
static char byte_from_wkb_state(wkb_parse_state *s)
{
    char char_value = 0;
    //LWDEBUG(4, "Entered function");

    wkb_parse_state_check(s, WKB_BYTE_SIZE);
    if (s->error)
        return 0;
    //LWDEBUG(4, "Passed state check");

    char_value = s->pos[0];
    //LWDEBUGF(4, "Read byte value: %x", char_value);
    s->pos += WKB_BYTE_SIZE;

    return char_value;
}


/**
* Int32
* Read 4-byte integer and advance the parse state forward.
*/
static uint32_t integer_from_wkb_state(wkb_parse_state *s)
{
    uint32_t i = 0;

    wkb_parse_state_check(s, WKB_INT_SIZE);
    if (s->error)
        return 0;

    memcpy(&i, s->pos, WKB_INT_SIZE);

    /* Swap? Copy into a stack-allocated integer. */
    if( s->swap_bytes )
    {
        int j = 0;
        uint8_t tmp;

        for( j = 0; j < WKB_INT_SIZE/2; j++ )
        {
            tmp = ((uint8_t*)(&i))[j];
            ((uint8_t*)(&i))[j] = ((uint8_t*)(&i))[WKB_INT_SIZE - j - 1];
            ((uint8_t*)(&i))[WKB_INT_SIZE - j - 1] = tmp;
        }
    }

    s->pos += WKB_INT_SIZE;
    return i;
}

/**
* Double
* Read an 8-byte double and advance the parse state forward.
*/
static double double_from_wkb_state(wkb_parse_state *s)
{
    double d = 0;

    memcpy(&d, s->pos, WKB_DOUBLE_SIZE);

    /* Swap? Copy into a stack-allocated integer. */
    if( s->swap_bytes )
    {
        int i = 0;
        uint8_t tmp;

        for( i = 0; i < WKB_DOUBLE_SIZE/2; i++ )
        {
            tmp = ((uint8_t*)(&d))[i];
            ((uint8_t*)(&d))[i] = ((uint8_t*)(&d))[WKB_DOUBLE_SIZE - i - 1];
            ((uint8_t*)(&d))[WKB_DOUBLE_SIZE - i - 1] = tmp;
        }

    }

    s->pos += WKB_DOUBLE_SIZE;
    return d;
}

/**
* Take in an unknown kind of wkb type number and ensure it comes out
* as an extended WKB type number (with Z/M/SRID flags masked onto the
* high bits).
*/
static void lwtype_from_wkb_state(wkb_parse_state *s, uint32_t wkb_type)
{
    uint32_t wkb_simple_type;

    s->has_z = LW_FALSE;
    s->has_m = LW_FALSE;
    s->has_srid = LW_FALSE;

    /* If any of the higher bits are set, this is probably an extended type. */
    if( wkb_type & 0xF0000000 )
    {
        if( wkb_type & WKBZOFFSET ) s->has_z = LW_TRUE;
        if( wkb_type & WKBMOFFSET ) s->has_m = LW_TRUE;
        if( wkb_type & WKBSRIDFLAG ) s->has_srid = LW_TRUE;
        LWDEBUGF(4, "Extended type: has_z=%d has_m=%d has_srid=%d", s->has_z, s->has_m, s->has_srid);
    }

    /* Mask off the flags */
    wkb_type = wkb_type & 0x0FFFFFFF;
    /* Strip out just the type number (1-12) from the ISO number (eg 3001-3012) */
    wkb_simple_type = wkb_type % 1000;

    /* Extract the Z/M information from ISO style numbers */
    if( wkb_type >= 3000 && wkb_type < 4000 )
    {
        s->has_z = LW_TRUE;
        s->has_m = LW_TRUE;
    }
    else if ( wkb_type >= 2000 && wkb_type < 3000 )
    {
        s->has_m = LW_TRUE;
    }
    else if ( wkb_type >= 1000 && wkb_type < 2000 )
    {
        s->has_z = LW_TRUE;
    }

    switch (wkb_simple_type)
    {
        case WKB_POINT_TYPE:
            s->lwtype = POINTTYPE;
            break;
        case WKB_LINESTRING_TYPE:
            s->lwtype = LINETYPE;
            break;
        case WKB_POLYGON_TYPE:
            s->lwtype = POLYGONTYPE;
            break;
        case WKB_MULTIPOINT_TYPE:
            s->lwtype = MULTIPOINTTYPE;
            break;
        case WKB_MULTILINESTRING_TYPE:
            s->lwtype = MULTILINETYPE;
            break;
        case WKB_MULTIPOLYGON_TYPE:
            s->lwtype = MULTIPOLYGONTYPE;
            break;
        case WKB_GEOMETRYCOLLECTION_TYPE:
            s->lwtype = COLLECTIONTYPE;
            break;
        case WKB_CIRCULARSTRING_TYPE:
            s->lwtype = CIRCSTRINGTYPE;
            break;
        case WKB_COMPOUNDCURVE_TYPE:
            s->lwtype = COMPOUNDTYPE;
            break;
        case WKB_CURVEPOLYGON_TYPE:
            s->lwtype = CURVEPOLYTYPE;
            break;
        case WKB_MULTICURVE_TYPE:
            s->lwtype = MULTICURVETYPE;
            break;
        case WKB_MULTISURFACE_TYPE:
            s->lwtype = MULTISURFACETYPE;
            break;
        case WKB_POLYHEDRALSURFACE_TYPE:
            s->lwtype = POLYHEDRALSURFACETYPE;
            break;
        case WKB_TIN_TYPE:
            s->lwtype = TINTYPE;
            break;
        case WKB_TRIANGLE_TYPE:
            s->lwtype = TRIANGLETYPE;
            break;

            /* PostGIS 1.5 emits 13, 14 for CurvePolygon, MultiCurve */
            /* These numbers aren't SQL/MM (numbers currently only */
            /* go up to 12. We can handle the old data here (for now??) */
            /* converting them into the lwtypes that are intended. */
        case WKB_CURVE_TYPE:
            s->lwtype = CURVEPOLYTYPE;
            break;
        case WKB_SURFACE_TYPE:
            s->lwtype = MULTICURVETYPE;
            break;

        default: /* Error! */
            //lwerror("Unknown WKB type (%d)! Full WKB type number was (%d).", wkb_simple_type, wkb_type);
            break;
    }

    return;
}

/**
* POINTARRAY
* Read a dynamically sized point array and advance the parse state forward.
* First read the number of points, then read the points.
*/
static POINTARRAY* ptarray_from_wkb_state(wkb_parse_state *s)
{
    POINTARRAY *pa = NULL;
    size_t pa_size;
    uint32_t ndims = 2;
    uint32_t npoints = 0;
    static uint32_t maxpoints = UINT_MAX / WKB_DOUBLE_SIZE / 4;

    /* Calculate the size of this point array. */
    npoints = integer_from_wkb_state(s);
    if (s->error)
        return NULL;
    if (npoints > maxpoints)
    {
        lwerror("Pointarray length (%d) is too large");
        return NULL;
    }

    if( s->has_z ) ndims++;
    if( s->has_m ) ndims++;
    pa_size = npoints * ndims * WKB_DOUBLE_SIZE;

    /* Empty! */
    if( npoints == 0 )
        return ptarray_construct(s->has_z, s->has_m, npoints);

    /* Does the data we want to read exist? */
    wkb_parse_state_check(s, pa_size);
    if (s->error)
        return NULL;

    /* If we're in a native endianness, we can just copy the data directly! */
    if( ! s->swap_bytes )
    {
        pa = ptarray_construct_copy_data(s->has_z, s->has_m, npoints, (uint8_t*)s->pos);
        s->pos += pa_size;
    }
        /* Otherwise we have to read each double, separately. */
    else
    {
        uint32_t i = 0;
        double *dlist;
        pa = ptarray_construct(s->has_z, s->has_m, npoints);
        dlist = (double*)(pa->serialized_pointlist);
        for( i = 0; i < npoints * ndims; i++ )
        {
            dlist[i] = double_from_wkb_state(s);
        }
    }

    return pa;
}

GBOX* gbox_copy(const GBOX *box)
{
    GBOX *copy = (GBOX*)lwalloc(sizeof(GBOX));
    memcpy(copy, box, sizeof(GBOX));
    return copy;
}

/**
 * @brief Deep clone a pointarray (also clones serialized pointlist)
 */
POINTARRAY *
ptarray_clone_deep(const POINTARRAY *in)
{
    POINTARRAY *out = lwalloc(sizeof(POINTARRAY));

    out->flags = in->flags;
    out->npoints = in->npoints;
    out->maxpoints = in->npoints;

    FLAGS_SET_READONLY(out->flags, 0);

    if (!in->npoints)
    {
        // Avoid calling lwalloc of 0 bytes
        out->serialized_pointlist = NULL;
    }
    else
    {
        size_t size = in->npoints * ptarray_point_size(in);
        out->serialized_pointlist = lwalloc(size);
        memcpy(out->serialized_pointlist, in->serialized_pointlist, size);
    }

    return out;
}

int
ptarray_is_closed_2d(const POINTARRAY *in)
{
    if (!in)
    {
        lwerror("ptarray_is_closed_2d: called with null point array");
        return 0;
    }
    if (in->npoints <= 1 ) return in->npoints; /* single-point are closed, empty not closed */

    return 0 == memcmp(getPoint_internal(in, 0), getPoint_internal(in, in->npoints-1), sizeof(POINT2D) );
}

/**
* Add a ring to a polygon. Point array will be referenced, not copied.
*/
int
lwpoly_add_ring(LWPOLY *poly, POINTARRAY *pa)
{
    if( ! poly || ! pa )
        return LW_FAILURE;

    /* We have used up our storage, add some more. */
    if( poly->nrings >= poly->maxrings )
    {
        int new_maxrings = 2 * (poly->nrings + 1);
        poly->rings = lwrealloc(poly->rings, new_maxrings * sizeof(POINTARRAY*));
        poly->maxrings = new_maxrings;
    }

    /* Add the new ring entry. */
    poly->rings[poly->nrings] = pa;
    poly->nrings++;

    return LW_SUCCESS;
}

/* Deep clone LWPOLY object. POINTARRAY are copied, as is ring array */
LWPOLY *
lwpoly_clone_deep(const LWPOLY *g)
{
    uint32_t i;
    LWPOLY *ret = lwalloc(sizeof(LWPOLY));
    memcpy(ret, g, sizeof(LWPOLY));
    if ( g->bbox ) ret->bbox = gbox_copy(g->bbox);
    ret->rings = lwalloc(sizeof(POINTARRAY *)*g->nrings);
    for ( i = 0; i < ret->nrings; i++ )
    {
        ret->rings[i] = ptarray_clone_deep(g->rings[i]);
    }
    FLAGS_SET_READONLY(ret->flags,0);
    return ret;
}

/* Deep clone LWLINE object. POINTARRAY *is* copied. */
LWLINE *
lwline_clone_deep(const LWLINE *g)
{
    LWLINE *ret = lwalloc(sizeof(LWLINE));

    LWDEBUGF(2, "lwline_clone_deep called with %p", g);
    memcpy(ret, g, sizeof(LWLINE));

    if ( g->bbox ) ret->bbox = gbox_copy(g->bbox);
    if ( g->points ) ret->points = ptarray_clone_deep(g->points);
    FLAGS_SET_READONLY(ret->flags,0);

    return ret;
}

/**
* @brief Deep clone #LWCOLLECTION object. #POINTARRAY are copied.
*/
LWCOLLECTION *
lwcollection_clone_deep(const LWCOLLECTION *g)
{
    uint32_t i;
    LWCOLLECTION *ret = lwalloc(sizeof(LWCOLLECTION));
    memcpy(ret, g, sizeof(LWCOLLECTION));
    if ( g->ngeoms > 0 )
    {
        ret->geoms = lwalloc(sizeof(LWGEOM *)*g->ngeoms);
        for (i=0; i<g->ngeoms; i++)
        {
            ret->geoms[i] = lwgeom_clone_deep(g->geoms[i]);
        }
        if ( g->bbox ) ret->bbox = gbox_copy(g->bbox);
    }
    else
    {
        ret->bbox = NULL; /* empty collection */
        ret->geoms = NULL;
    }
    return ret;
}

/**
* Deep-clone an #LWGEOM object. #POINTARRAY <em>are</em> copied.
*/
LWGEOM *
lwgeom_clone_deep(const LWGEOM *lwgeom)
{

    switch (lwgeom->type)
    {
        case POINTTYPE:
        case LINETYPE:
        case CIRCSTRINGTYPE:
        case TRIANGLETYPE:
            return (LWGEOM *)lwline_clone_deep((LWLINE *)lwgeom);
        case POLYGONTYPE:
            return (LWGEOM *)lwpoly_clone_deep((LWPOLY *)lwgeom);
        case COMPOUNDTYPE:
        case CURVEPOLYTYPE:
        case MULTICURVETYPE:
        case MULTISURFACETYPE:
        case MULTIPOINTTYPE:
        case MULTILINETYPE:
        case MULTIPOLYGONTYPE:
        case POLYHEDRALSURFACETYPE:
        case TINTYPE:
        case COLLECTIONTYPE:
            return (LWGEOM *)lwcollection_clone_deep((LWCOLLECTION *)lwgeom);
        default:
            lwerror("lwgeom_clone_deep: Unknown geometry type: %s", lwtype_name(lwgeom->type));
            return NULL;
    }
}

/**
* POINT
* Read a WKB point, starting just after the endian byte,
* type number and optional srid number.
* Advance the parse state forward appropriately.
* WKB point has just a set of doubles, with the quantity depending on the
* dimension of the point, so this looks like a special case of the above
* with only one point.
*/
static LWPOINT* lwpoint_from_wkb_state(wkb_parse_state *s)
{
    static uint32_t npoints = 1;
    POINTARRAY *pa = NULL;
    size_t pa_size;
    uint32_t ndims = 2;
    const POINT2D *pt;

    /* Count the dimensions. */
    if( s->has_z ) ndims++;
    if( s->has_m ) ndims++;
    pa_size = ndims * WKB_DOUBLE_SIZE;

    /* Does the data we want to read exist? */
    wkb_parse_state_check(s, pa_size);
    if (s->error)
        return NULL;

    /* If we're in a native endianness, we can just copy the data directly! */
    if( ! s->swap_bytes )
    {
        pa = ptarray_construct_copy_data(s->has_z, s->has_m, npoints, (uint8_t*)s->pos);
        s->pos += pa_size;
    }
        /* Otherwise we have to read each double, separately */
    else
    {
        uint32_t i = 0;
        double *dlist;
        pa = ptarray_construct(s->has_z, s->has_m, npoints);
        dlist = (double*)(pa->serialized_pointlist);
        for( i = 0; i < ndims; i++ )
        {
            dlist[i] = double_from_wkb_state(s);
        }
    }

    /* Check for POINT(NaN NaN) ==> POINT EMPTY */
    pt = getPoint2d_cp(pa, 0);
    if ( isnan(pt->x) && isnan(pt->y) )
    {
        ptarray_free(pa);
        return lwpoint_construct_empty(s->srid, s->has_z, s->has_m);
    }
    else
    {
        return lwpoint_construct(s->srid, NULL, pa);
    }
}

/**
* POLYGON
* Read a WKB polygon, starting just after the endian byte,
* type number and optional srid number. Advance the parse state
* forward appropriately.
* First read the number of rings, then read each ring
* (which are structured as point arrays)
*/
static LWPOLY* lwpoly_from_wkb_state(wkb_parse_state *s)
{
    LWPOLY *poly;
    uint32_t i = 0;
    uint32_t nrings = integer_from_wkb_state(s);
    if (s->error) {
        return NULL;
    }

    poly = lwpoly_construct_empty(s->srid, s->has_z, s->has_m);

    LWDEBUGF(4,"Polygon has %d rings", nrings);

    /* Empty polygon? */
    if( nrings == 0 )
        return poly;

    for( i = 0; i < nrings; i++ )
    {
        POINTARRAY *pa = ptarray_from_wkb_state(s);
        if (pa == NULL)
        {
            lwpoly_free(poly);
            return NULL;
        }

        /* Check for at least four points. */
        if (s->check & LW_PARSER_CHECK_MINPOINTS && pa->npoints < 4)
        {
            lwpoly_free(poly);
            LWDEBUGF(2, "%s must have at least four points in each ring", lwtype_name(s->lwtype));
            lwerror("%s must have at least four points in each ring", lwtype_name(s->lwtype));
            return NULL;
        }

        /* Check that first and last points are the same. */
        if( s->check & LW_PARSER_CHECK_CLOSURE && ! ptarray_is_closed_2d(pa) )
        {
            lwpoly_free(poly);
            LWDEBUGF(2, "%s must have closed rings", lwtype_name(s->lwtype));
            lwerror("%s must have closed rings", lwtype_name(s->lwtype));
            return NULL;
        }

        /* Add ring to polygon */
        if ( lwpoly_add_ring(poly, pa) == LW_FAILURE )
        {
            lwpoly_free(poly);
            LWDEBUG(2, "Unable to add ring to polygon");
            lwerror("Unable to add ring to polygon");
            return NULL;
        }

    }
    return poly;
}

/**
* CURVEPOLYTYPE
*/
static LWPOLY* lwcurvepoly_from_wkb_state(wkb_parse_state *s)
{
    LWCURVEPOLY *cp;
    LWGEOM *geom;
    uint32_t i;
    uint32_t ngeoms = integer_from_wkb_state(s);
    if (s->error)
        return NULL;
    cp = lwcurvepoly_construct_empty(s->srid, s->has_z, s->has_m);
    geom = NULL;


    /* Empty collection? */
    if ( ngeoms == 0 )
        return (LWPOLY*)cp;

    for ( i = 0; i < ngeoms; i++ )
    {
        geom = lwgeom_from_wkb_state(s);
        if ( lwcurvepoly_add_ring(cp, geom) == LW_FAILURE )
        {
            lwgeom_free(geom);
            lwgeom_free((LWGEOM *)cp);
            lwerror("Unable to add geometry (%p) to curvepoly (%p)", geom, cp);
            return NULL;
        }
    }

    return (LWPOLY*)cp;
}

/**
* Appends geom to the collection managed by col. Does not copy or
* clone, simply takes a reference on the passed geom.
*/
LWCOLLECTION* lwcollection_add_lwgeom(LWCOLLECTION *col, const LWGEOM *geom)
{
    if (!col || !geom) return NULL;

    if (!col->geoms && (col->ngeoms || col->maxgeoms))
    {
        lwerror("Collection is in inconsistent state. Null memory but non-zero collection counts.");
        return NULL;
    }

    /* Check type compatibility */
    if ( ! lwcollection_allows_subtype(col->type, geom->type) ) {
        lwerror("%s cannot contain %s element", lwtype_name(col->type), lwtype_name(geom->type));
        return NULL;
    }

    /* In case this is a truly empty, make some initial space  */
    if (!col->geoms)
    {
        col->maxgeoms = 2;
        col->ngeoms = 0;
        col->geoms = lwalloc(col->maxgeoms * sizeof(LWGEOM*));
    }

    /* Allocate more space if we need it */
    lwcollection_reserve(col, col->ngeoms + 1);

#if PARANOIA_LEVEL > 1
    /* See http://trac.osgeo.org/postgis/ticket/2933 */
	/* Make sure we don't already have a reference to this geom */
	{
		uint32_t i = 0;
		for (i = 0; i < col->ngeoms; i++)
		{
			if (col->geoms[i] == geom)
			{
				lwerror("%s [%d] found duplicate geometry in collection %p == %p", __FILE__, __LINE__, col->geoms[i], geom);
				return col;
			}
		}
	}
#endif

    col->geoms[col->ngeoms] = (LWGEOM*)geom;
    col->ngeoms++;
    return col;
}

/**
* POLYHEDRALSURFACETYPE
*/

/**
* COLLECTION, MULTIPOINTTYPE, MULTILINETYPE, MULTIPOLYGONTYPE, COMPOUNDTYPE,
* MULTICURVETYPE, MULTISURFACETYPE,
* TINTYPE
*/
static LWCOLLECTION* lwcollection_from_wkb_state(wkb_parse_state *s)
{
    LWCOLLECTION *col;
    LWGEOM *geom;
    uint32_t i;
    uint32_t ngeoms = integer_from_wkb_state(s);
    if (s->error)
        return NULL;

    col = lwcollection_construct_empty(s->lwtype, s->srid, s->has_z, s->has_m);
    geom = NULL;


    /* Empty collection? */
    if ( ngeoms == 0 )
        return col;

    /* Be strict in polyhedral surface closures */
    //if ( s->lwtype == POLYHEDRALSURFACETYPE )
    //    s->check |= LW_PARSER_CHECK_ZCLOSURE;

    for ( i = 0; i < ngeoms; i++ )
    {
        geom = lwgeom_from_wkb_state(s);
        if ( lwcollection_add_lwgeom(col, geom) == NULL )
        {
            lwgeom_free(geom);
            lwgeom_free((LWGEOM *)col);
            lwerror("Unable to add geometry (%p) to collection (%p)", geom, col);
            return NULL;
        }
    }

    return col;
}

/**
* GEOMETRY
* Generic handling for WKB geometries. The front of every WKB geometry
* (including those embedded in collections) is an endian byte, a type
* number and an optional srid number. We handle all those here, then pass
* to the appropriate handler for the specific type.
*/
LWGEOM* lwgeom_from_wkb_state(wkb_parse_state *s)
{
    char wkb_little_endian;
    uint32_t wkb_type;

    /* Fail when handed incorrect starting byte */
    wkb_little_endian = byte_from_wkb_state(s);
    if (s->error)
        return NULL;
    if( wkb_little_endian != 1 && wkb_little_endian != 0 )
    {
        LWDEBUG(4,"Leaving due to bad first byte!");
        lwerror("Invalid endian flag value encountered. =%c", wkb_little_endian);
        return NULL;
    }

    /* Check the endianness of our input  */
    s->swap_bytes = LW_FALSE;

    /* Machine arch is big endian, request is for little */
    if (IS_BIG_ENDIAN && wkb_little_endian)
        s->swap_bytes = LW_TRUE;
        /* Machine arch is little endian, request is for big */
    else if ((!IS_BIG_ENDIAN) && (!wkb_little_endian))
        s->swap_bytes = LW_TRUE;

    /* Read the type number */
    wkb_type = integer_from_wkb_state(s);
    if (s->error)
        return NULL;
    lwtype_from_wkb_state(s, wkb_type);

    /* Read the SRID, if necessary */
    if( s->has_srid )
    {
        s->srid = clamp_srid(integer_from_wkb_state(s));
        if (s->error)
            return NULL;
        /* TODO: warn on explicit UNKNOWN srid ? */
        LWDEBUGF(4,"Got SRID: %u", s->srid);
    }

    /* Do the right thing */
    switch( s->lwtype )
    {
        case POINTTYPE:
            return (LWGEOM*)lwpoint_from_wkb_state(s);
            break;

        case POLYGONTYPE:
            return (LWGEOM*)lwpoly_from_wkb_state(s);
            break;

        case CURVEPOLYTYPE:
            return (LWGEOM*)lwcurvepoly_from_wkb_state(s);
            break;
        case MULTIPOINTTYPE:
        case MULTIPOLYGONTYPE:
        case COLLECTIONTYPE:
            return (LWGEOM*)lwcollection_from_wkb_state(s);
            break;
        case LINETYPE:
        case CIRCSTRINGTYPE:
        case TRIANGLETYPE:
        case MULTICURVETYPE:
        case MULTISURFACETYPE:
        case COMPOUNDTYPE:
        case MULTILINETYPE:
        case POLYHEDRALSURFACETYPE:
        case TINTYPE:
        default:
            lwerror("Unsupported geometry type: %s", lwtype_name(s->lwtype));
    }

    /* Return value to keep compiler happy. */
    return NULL;

}

/*
 * Copies a point from the point array into the parameter point
 * will set point's z=NO_Z_VALUE if pa is 2d
 * will set point's m=NO_M_VALUE if pa is 3d or 2d
 *
 * NOTE: point is a real POINT3D *not* a pointer
 */
POINT4D
getPoint4d(const POINTARRAY *pa, uint32_t n)
{
    POINT4D result;
    getPoint4d_p(pa, n, &result);
    return result;
}

/*
 * Copies a point from the point array into the parameter point
 * will set point's z=NO_Z_VALUE  if pa is 2d
 * will set point's m=NO_M_VALUE  if pa is 3d or 2d
 *
 * NOTE: this will modify the point4d pointed to by 'point'.
 *
 * @return 0 on error, 1 on success
 */
int
getPoint4d_p(const POINTARRAY *pa, uint32_t n, POINT4D *op)
{
    uint8_t *ptr;
    int zmflag;

    if ( ! pa )
    {
        lwerror("%s [%d] NULL POINTARRAY input", __FILE__, __LINE__);
        return 0;
    }

    if ( n>=pa->npoints )
    {
        lwnotice("%s [%d] called with n=%d and npoints=%d", __FILE__, __LINE__, n, pa->npoints);
        return 0;
    }

    LWDEBUG(5, "getPoint4d_p called.");

    /* Get a pointer to nth point offset and zmflag */
    ptr=getPoint_internal(pa, n);
    zmflag=FLAGS_GET_ZM(pa->flags);

    LWDEBUGF(5, "ptr %p, zmflag %d", ptr, zmflag);

    switch (zmflag)
    {
        case 0: /* 2d  */
            memcpy(op, ptr, sizeof(POINT2D));
            op->m=NO_M_VALUE;
            op->z=NO_Z_VALUE;
            break;

        case 3: /* ZM */
            memcpy(op, ptr, sizeof(POINT4D));
            break;

        case 2: /* Z */
            memcpy(op, ptr, sizeof(POINT3DZ));
            op->m=NO_M_VALUE;
            break;

        case 1: /* M */
            memcpy(op, ptr, sizeof(POINT3DM));
            op->m=op->z; /* we use Z as temporary storage */
            op->z=NO_Z_VALUE;
            break;

        default:
            lwerror("Unknown ZM flag ??");
            return 0;
    }
    return 1;

}
