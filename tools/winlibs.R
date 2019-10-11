if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# For details see: https://github.com/rwinlib/gdal2
VERSION <- commandArgs(TRUE)
if(!file.exists(sprintf("../windows/gdal2-%s/include/gdal/gdal.h", VERSION))){
  if(getRversion() < "3.3.0") setInternet2()
  dir.create("../windows", showWarnings = FALSE)
  download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = TRUE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
  download.file(sprintf("https://github.com/rgeoda/libgeoda_static/archive/0.0.3.3.zip", VERSION), "geodalib.zip", quiet = TRUE)
  unzip("geodalib.zip", exdir = "../windows")
  unlink("geodalib.zip")
  download.file(sprintf("https://codeload.github.com/rgeoda/wx_static/zip/master", VERSION), "wxlib.zip", quiet = TRUE)
  unzip("wxlib.zip", exdir = "../windows")
  unlink("wxlib.zip")
  download.file(sprintf("https://github.com/rgeoda/ANN_static/archive/v1.1.2.zip", VERSION), "annlib.zip", quiet = TRUE)
  unzip("annlib.zip", exdir = "../windows")
  unlink("annlib.zip")
  download.file(sprintf("https://github.com/rgeoda/boost_static/archive/v1.57.0.zip", VERSION), "boostlib.zip", quiet = TRUE)
  unzip("boostlib.zip", exdir = "../windows")
  unlink("boostlib.zip")
}
