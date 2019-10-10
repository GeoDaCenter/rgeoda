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
  download.file("https://github.com/rgeoda/libgeoda_static/archive/0.0.3.2.zip", "geodalib.zip", quiet = TRUE)
  unzip("geodalib.zip", exdir = "../deps")
  unlink("geodalib.zip")
  download.file("https://codeload.github.com/rgeoda/wx_static/zip/master", "wxlib.zip", quiet = TRUE)
  unzip("wxlib.zip", exdir = "../deps")
  unlink("wxlib.zip")
  download.file("https://codeload.github.com/rgeoda/ANN_static/zip/master", "annlib.zip", quiet = TRUE)
  unzip("annlib.zip", exdir = "../deps")
  unlink("annlib.zip")
  download.file("https://codeload.github.com/rgeoda/boost_static/zip/master", "boostlib.zip", quiet = TRUE)
  unzip("boostlib.zip", exdir = "../deps")
  unlink("boostlib.zip")
}
