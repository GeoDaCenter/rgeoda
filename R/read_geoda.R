# create a nice wrapper R class for libgeoda-r
rgeoda <- setRefClass("rgeoda",
  fields = list(
    is_memory = "logical"
  ),
  slots=list(
    data="data.frame",
    geometries="list",
    geoda="_p_GeoDa"
  ),
  methods = list(
    CreateQueenWeights = function(...) {
    },
    CreateRookWeights = function(...) {
    },
    CreateKnnWeights = function(...) {
    },
    CreateKernelWeights = function(...) {
    },
    CreateSocialWeights = function(...) {
    }
  )
)

create_geoda_from_sf = function(obj) {
  print ("create geoda from sf object")
  
}

create_geoda_from_sp = function(obj) {
  print ("create geoda from sf object")

}

read_geoda = function(obj) {
  if (is(obj, "sf")) {
    create_geoda_from_sf(obj)

  } else if (is(obj, "sp")) {
    create_geoda_from_sp(obj)

  } else if (typeof(obj) == "character") {
    return <- GeoDa(obj)
  } else {
    stop ("Only a sf object or a string of input datasource is allowed")
  }
}