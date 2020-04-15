# get Island Scrub Jay data

# Data were collected at 307 survey locations ("point transects") on Santa Cruz Island, California during the Fall of 2008. The distance data are binned into 3 distance intervals [0-100], (100-200], and (200-300]. The coordinates of the survey locations as well as 3 habitat covariates are also included. 
library(unmarked)
data(issj)


# spring data here:
# https://figshare.com/articles/Supplement_1_R_code_data_and_grid_covariates_used_in_the_analyses_/3517754

# setup a Sample.Label
issj$Sample.Label <- 1:nrow(issj)

# standardize covars as in paper
s_elevation <- scale(issj$elevation)
s_chaparral <- scale(issj$chaparral)
s_forest <- scale(issj$forest)
s_elevation_detail <- c(attr(s_elevation, "scaled:center"),
                        attr(s_elevation, "scaled:scale"))
s_chaparral_detail <- c(attr(s_chaparral, "scaled:center"),
                        attr(s_chaparral, "scaled:scale"))
s_forest_detail <- c(attr(s_forest, "scaled:center"),
                        attr(s_forest, "scaled:scale"))
issj$elevation <- s_elevation[,1]
issj$chaparral <- s_chaparral[,1]
issj$forest <- s_forest[,1]


# get segment data
segs <- issj[, c("x", "y", "elevation", "forest",
                 "chaparral", "Sample.Label")]

# get observations
obs <- list()

mids <- c(50, 150, 250)
begs <- c(0, 100, 200)
ends <- c(100, 200, 300)

for(i in 1:nrow(issj)){
  this_row <- issj[i, ]

  if(all(this_row[1:3] == 0)) next


  obs[[i]] <- data.frame(distance  = rep(mids, this_row[1:3]),
                         distbegin = rep(begs, this_row[1:3]),
                         distend   = rep(ends, this_row[1:3]))
  obs[[i]]$forest    <- this_row$forest
  obs[[i]]$chaparral <- this_row$chaparral
  obs[[i]]$habitat <- this_row$habitat
  obs[[i]]$Sample.Label <- this_row$Sample.Label

}
# mudge back to a data.frame
obs <- do.call(rbind.data.frame, obs)

# add extra cols
obs$object <- 1:nrow(obs)
obs$size <- 1
segs$Effort <- 1

# that's the fall data done!
obs_fall <- obs

# get the spring data
source("data_100m.R")

obs <- list()
for(i in 1:nrow(Xall.spring.100m)){
  this_row <- Xall.spring.100m[i, ]

  if(all(this_row == 0)) next


  obs[[i]] <- data.frame(distance  = rep(mids, this_row),
                         distbegin = rep(begs, this_row),
                         distend   = rep(ends, this_row))
  obs[[i]]$forest    <- issj$forest[i]
  obs[[i]]$chaparral <- issj$chaparral[i]
  obs[[i]]$habitat <- issj$habitat[i]
  obs[[i]]$Sample.Label <- issj$Sample.Label[i]
}
# mudge back to a data.frame
obs <- do.call(rbind.data.frame, obs)

# add extra cols
obs$object <- 1:nrow(obs)
obs$size <- 1

# that's the fall data done!
obs_spring <- obs

# prediction grid
source("grid_covariates.R")
pred <- data.frame(chaparral = gchap,
                   elevation = gelev,
                   forest    = gfor,
                   off.set   = 300*300)
data(cruz)
cruz$off.set <- 300*300

save(obs_fall, obs_spring, segs, pred, cruz,
     s_chaparral_detail, s_forest_detail, s_elevation_detail,
     file="issj.RData")



