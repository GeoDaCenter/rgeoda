FROM rocker/r-base

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y git libssl-dev libgeos-dev libgeos++-dev gdal-bin libproj-dev libgdal-dev libudunits2-dev

RUN install2.r --error proxy Rcpp wk sp digest sf BH wkb TinyTex
