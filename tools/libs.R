if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# For details see: https://github.com/rwinlib/gdal2
VERSION <- commandArgs(TRUE)
if(!file.exists("../deps/gdal_static-master/include/gdal.h")){
  if(getRversion() < "3.3.0") setInternet2()
  dir.create("../deps", showWarnings = FALSE)
  download.file("https://codeload.github.com/libgeoda/gdal_static/zip/master", "lib.zip", quiet = TRUE)
  unzip("lib.zip", exdir = "../deps")
  unlink("lib.zip")
  download.file("https://github.com/rgeoda/libgeoda_static/releases/download/0.0.3.3/0.0.3.4.zip", "geodalib.zip", quiet = TRUE)
  unzip("geodalib.zip", exdir = "../deps")
  unlink("geodalib.zip")
  download.file("https://codeload.github.com/rgeoda/boost_static/zip/master", "boostlib.zip", quiet = TRUE)
  unzip("boostlib.zip", exdir = "../deps")
  unlink("boostlib.zip")
}
