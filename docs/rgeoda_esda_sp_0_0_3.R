## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
# use the Guerry.shp comes with the rgeoda package
guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")

# load rgdal library which sp is working with
library(rgdal)

## ------------------------------------------------------------------------
guerry_sp <- readOGR(guerry_path)

## ------------------------------------------------------------------------
plot(guerry_sp)

## ------------------------------------------------------------------------
library(rgeoda)

guerry <- sp_to_geoda(guerry_sp)

## ------------------------------------------------------------------------
queen_w <- queen_weights(guerry)
crm_prp <- as.numeric(as.character(guerry_sp$Crm_prp))

## ------------------------------------------------------------------------
lisa <- local_moran(queen_w, crm_prp)

## ------------------------------------------------------------------------
lisa_colors <- lisa$GetColors() 
lisa_labels <- lisa$GetLabels()
lisa_clusters <- lisa$GetClusterIndicators()

plot(guerry_sp, 
     col=sapply(lisa_clusters, function(x){return(lisa_colors[[x+1]])}), 
     border = "#333333", lwd=0.2)
title(main = "Local Moran Map of Crm_prp")
legend('bottomleft', legend = lisa_labels, fill = lisa_colors, border = "#eeeeee")

## ------------------------------------------------------------------------
guerry_sp$moran_cluster <- lisa_clusters

## ------------------------------------------------------------------------
lisa_clusters

## ------------------------------------------------------------------------
lisa_p <- lisa$GetPValues()
p_labels <- c("Not significant", "p <= 0.05", "p <= 0.01", "p <= 0.001")
p_colors <- c("#eeeeee", "#84f576", "#53c53c", "#348124")
plot(guerry_sp, 
     col=sapply(lisa_p, function(x){
       if (x <= 0.001) return(p_colors[4])
       else if (x <= 0.01) return(p_colors[3])
       else if (x <= 0.05) return (p_colors[2])
       else return(p_colors[1])
       }), 
     border = "#333333", lwd=0.2)
title(main = "Local Moran Map of Crm_prp")
legend('bottomleft', legend = p_labels, fill = p_colors, border = "#eeeeee")

## ------------------------------------------------------------------------
lisa <- local_geary(queen_w, crm_prp)
lisa_colors <- lisa$GetColors() 
lisa_labels <- lisa$GetLabels()
lisa_clusters <- lisa$GetClusterIndicators()

plot(guerry_sp, 
     col=sapply(lisa_clusters, function(x){return(lisa_colors[[x+1]])}), 
     border = "#333333", lwd=0.2)
title(main = "Local Geary Map of Crm_prp")
legend('bottomleft', legend = lisa_labels, fill = lisa_colors, border = "#eeeeee")

## ------------------------------------------------------------------------
lisa <- local_g(queen_w, crm_prp)
lisa_colors <- lisa$GetColors() 
lisa_labels <- lisa$GetLabels()
lisa_clusters <- lisa$GetClusterIndicators()

plot(guerry_sp, 
     col=sapply(lisa_clusters, function(x){return(lisa_colors[[x+1]])}), 
     border = "#333333", lwd=0.2)
title(main = "Local Getis-Ord's G Map of Crm_prp")
legend('bottomleft', legend = lisa_labels, fill = lisa_colors, border = "#eeeeee")

## ------------------------------------------------------------------------
Crm_prp <- as.numeric(as.character(guerry_sp$Crm_prp))
Litercy <- as.numeric(as.character(guerry_sp$Litercy))
Donatns <- as.numeric(as.character(guerry_sp$Donatns))
Infants <- as.numeric(as.character(guerry_sp$Infants))
Suicids <- as.numeric(as.character(guerry_sp$Suicids))
data <- list(crm_prp, Crm_prp, Donatns, Infants, Suicids)
guerry_clusters <- skater(4, queen_w, data)

# Get some colors for each clusters
skater_colors <- palette()[2:5]
skater_labels <- c("c1","c2","c3","c4")

# Assign a color for each observation
colors <- rep("#000000", queen_w$num_obs)
for (i in 1:4) {
  for (j in guerry_clusters[i]) {
    colors[j+1] <- skater_colors[i]
  }
}

# plot
plot(guerry_sp,  col=colors, border = "#333333", lwd=0.2)
title(main = "SKATER Clustering Map")
legend('bottomleft', legend = skater_labels, fill = skater_colors, border = "#eeeeee")

## ------------------------------------------------------------------------
bound_vals <- as.numeric(guerry_sp$Pop1831)
min_bound <- 3236.67 # 10% of Pop1831

maxp_clusters <- maxp(queen_w, data, bound_vals, min_bound, "greedy")

# Get some colors for each clusters
maxp_colors <- palette()[2:10]
maxp_labels <- c("c1","c2","c3","c4","c5","c6","c7","c8")

# Assign a color for each observation
colors <- rep("#000000", queen_w$num_obs)
for (i in 1:8) {
  for (j in maxp_clusters[i]) {
    colors[j+1] <- maxp_colors[i]
  }
}

# plot
plot(guerry_sp,  col=colors, border = "#333333", lwd=0.2)
title(main = "Max-p Clustering Map")
legend('bottomleft', legend = maxp_labels, fill = maxp_colors, border = "#eeeeee")

