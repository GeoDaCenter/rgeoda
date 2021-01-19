#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h> /* for tolower */

#include "config.h"
#include "utils.h"


/* Default reporters */
static void default_noticereporter(const char *fmt, va_list ap);
static void default_errorreporter(const char *fmt, va_list ap);
//lwreporter lwnotice_var = default_noticereporter;
//lwreporter lwerror_var = default_errorreporter;

/* Default logger */
static void default_debuglogger(int level, const char *fmt, va_list ap);
//lwdebuglogger lwdebug_var = default_debuglogger;

void
lwnotice(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);

    /* Call the supplied function */
    default_noticereporter(fmt, ap);

    va_end(ap);
}

void
lwerror(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);

    /* Call the supplied function */
   default_errorreporter(fmt, ap);

    va_end(ap);
}

void
lwdebug(int level, const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);

    /* Call the supplied function */
    default_debuglogger(level, fmt, ap);

    va_end(ap);
}

static void
default_debuglogger(int level, const char *fmt, va_list ap)
{
#ifdef __PGGEODA__
    char msg[LW_MSG_MAXLEN+1];
    if ( POSTGIS_DEBUG_LEVEL >= level )
    {
        /* Space pad the debug output */
        int i;
        for ( i = 0; i < level; i++ )
            msg[i] = ' ';
        vsnprintf(msg+i, LW_MSG_MAXLEN-i, fmt, ap);
        msg[LW_MSG_MAXLEN]='\0';
        fprintf(stderr, "%s\n", msg);
    }
#endif
}

static void
default_errorreporter(const char *fmt, va_list ap)
{
#ifdef __PGGEODA__
    char msg[LW_MSG_MAXLEN+1];
    vsnprintf (msg, LW_MSG_MAXLEN, fmt, ap);
    msg[LW_MSG_MAXLEN]='\0';
    fprintf(stderr, "%s\n", msg);
    exit(1);
#endif
}


static void
default_noticereporter(const char *fmt, va_list ap)
{
#ifdef __PGGEODA__
    char msg[LW_MSG_MAXLEN+1];
    vsnprintf (msg, LW_MSG_MAXLEN, fmt, ap);
    msg[LW_MSG_MAXLEN]='\0';
    fprintf(stderr, "%s\n", msg);
#endif
}
