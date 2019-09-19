AC_INIT([rgeodap], 0.0.1)

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

BUILD_PLATFORM=unknown
WINDOWS_BUILD=no
OS_NAME=linux
AC_MSG_CHECKING([for build platform])
# see https://stackoverflow.com/questions/714100/os-detecting-makefile
# but we do this in the configure for homogeneity
case "${OS}" in
    Windows_NT*)
        BUILD_PLATFORM=Windows
        ;;
    *)
        uname_str=`uname -s`
        case "${uname_str}" in
            Linux*)
                BUILD_PLATFORM=Linux
                OS_NAME=linux
                ;;
            Darwin*)
                BUILD_PLATFORM=Darwin
                OS_NAME=osx
                ;;
            cygwin*|mingw*)
                BUILD_PLATFORM=Windows
                OS_NAME=mingw64
                ;;
            *)
                AC_MSG_ERROR(["OS $unamestr is not supported"])
                ;;
        esac
esac
AC_MSG_RESULT([$BUILD_PLATFORM])
AC_SUBST([BUILD_PLATFORM])
AC_SUBST([OS_NAME])

# default values
AC_SUBST([GEODA_INCLUDE], "-I./deps/libgeoda/include -I./deps/include -I./deps/lib/linux/wx/include/base-unicode-static-3.0 -D_FILE_OFFSET_BITS=64")
AC_SUBST([GEODA_LIBS], "./deps/libgeoda/lib/linux/libgeoda.a ./deps/lib/linux/libgdal.a ./deps/lib/linux/libgeos_c.a ./deps/lib/linux/libgeos.a ./deps/lib/linux/libANN.a  ./deps/lib/linux/libwx_baseu-3.0.a ./deps/lib/linux/libwxregexu-3.0.a ./deps/lib/linux/libboost_system-mt.a ./deps/lib/linux/libboost_chrono-mt.a ./deps/lib/linux/libboost_date_time-mt.a ./deps/lib/linux/libboost_thread-mt.a -L/usr/lib -dt ./deps/lib/linux/libiconv.a")

AS_IF([test x"$BUILD_PLATFORM" = x"Darwin"], [
    GEODA_INCLUDE="-I./deps/libgeoda/include -I./deps/include -I./deps/lib/osx/wx/include/base-unicode-static-3.0 -D_FILE_OFFSET_BITS=64"
	GEODA_LIBS="-L/usr/lib -lcurl -liconv ./deps/libgeoda/lib/osx/libgeoda.a ./deps/lib/osx/libANN.a ./deps/lib/osx/libgdal.a ./deps/lib/osx/libwx_baseu-3.0.a ./deps/lib/osx/libwxregexu-3.0.a ./deps/lib/osx/libgeos.a ./deps/lib/osx/libgeos_c.a -framework IOKit -framework CoreServices -framework System -framework ApplicationServices ./deps/lib/osx/libboost_thread.a ./deps/lib/osx/libboost_system.a ./deps/lib/osx/libboost_chrono.a ./deps/lib/osx/libboost_date_time.a"
])

# allow for override
AC_ARG_WITH([geoda],
            AC_HELP_STRING([--with-geoda=PREFIX],
                           [path to where geoda is installed]),
            [
             GEODA_INCLUDE="-I${with_geoda}/include"
             GEODA_LIBS="${with_geoda}/lib/libgeoda.a"
            ],
            [])

# create and report output
AS_IF([test x"$BUILD_PLATFORM" != x"Windows"], [
    AC_CONFIG_COMMANDS_PRE([mv _src src])
])
AS_IF([test x"$BUILD_PLATFORM" != x"Windows"], [
    AC_CONFIG_FILES([src/Makevars])
    AC_OUTPUT
])