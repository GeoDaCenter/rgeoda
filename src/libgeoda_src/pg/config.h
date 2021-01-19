#ifndef POSTGEODA_CONFIG_H
#define POSTGEODA_CONFIG_H

#ifndef DEBUG
/* #undef DEBUG */
#endif

#define POSTGIS_DEBUG_LEVEL 1

/* Define to 1 to enable memory checks in pointarray management. */
#undef PARANOIA_LEVEL

/**
* Maximum allowed SRID value in serialized geometry.
* Currently we are using 21 bits (2097152) of storage for SRID.
*/
#define SRID_MAXIMUM 999999

/**
 * Maximum valid SRID value for the user
 * We reserve 1000 values for internal use
 */
#define SRID_USER_MAXIMUM 998999

/** Unknown SRID value */
#define SRID_UNKNOWN 0
#define SRID_IS_UNKNOWN(x) ((int)x<=0)

#endif
