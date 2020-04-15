# school size with varprop
# factor-smooth interaction with chopped-up school sizes
# example using harbour porpoises from parts of SCANS-II

# libraries and helper functions
library(mrds)
library(dsm)
library(plyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("R/make_ss_fs_data.R")
source("R/est_N_var_ss.R")


## 1. data setup

# load the SCANS-II data
load("data/SCANSII.RData")

# select the area we are interested in:
#   R is coastal western Ireland
#   O is the Irish sea
#   N is Western coastal Scotland
area_of_interest <- c("N", "R", "O")
air.ddf <- air.ddf[air.ddf$Stratum %in% area_of_interest, ]
air.segments <- air.segments[air.segments$Stratum %in% area_of_interest, ]

# select the harbour porpoise observations only
air.ddf <- air.ddf[air.ddf$Species=="ppho", ]

# only include segments with beaufort <= max.bf
max.bf <- 4
air.segments <- air.segments[air.segments$Beaufort <= max.bf, ]
air.ddf <- air.ddf[air.ddf$Beaufort <= max.bf, ]

# shortcut names
ddf_data <- air.ddf
segments_data <- air.segments

# bin Beaufort
ddf_data$beaufort <- cut(ddf_data$Beaufort, c(0,1,2,4), include.lowest=TRUE)
segments_data$beaufort <- cut(segments_data$Beaufort, c(0,1,2,4),
                              include.lowest=TRUE)



## 2. detection function fitting, looking at group size relationship
width <- 300
df_test <- ddf(method='ds',
               dsmodel=~mcds(formula=~size+beaufort, key="hr"),
               data=ddf_data, meta.data=list(width=width))
ddf_data <- ddf_data[ddf_data$distance <= width, ]

plot(ddf_data$size, predict(df_test)[[1]], xlab="size", ylab="p")


## 3. chop-up the data
cut.labs <- c("Singletons", "Twos", "3-5")
fs_data <- make_ss_fs_data(ddf_data, segments_data,
                           c(0, 1, 2, 5),
                           cut.labs)
ddf_data_ss <- fs_data$obs
segments_data_ss <- fs_data$segs



## 4. fit detection function to new factor school size
df_s <- ddf(method='ds',
            dsmodel=~mcds(formula=~size_class+beaufort, key="hr"),
            data=ddf_data_ss, meta.data=list(width=width))



## 5. fit some spatial models

# dsm with school size-xy "interaction"
df_fs_chop <- dsm(count ~ s(x, y, size_class, k=20, bs="fs"),
                  ddf.obj = df_s, segment.data = segments_data_ss,
                  observation.data = ddf_data_ss, family=tw(a=1.2))

# how well does the model do?
source("R/obs_exp.R")
print(obs_exp(df_fs_chop, ~beaufort))


## comparison with estimated abundance model
dsm_boring <- dsm(abundance.est ~ s(x, y, k=60),
          ddf.obj = df_s, segment.data = segments_data,
          observation.data = ddf_data, family=tw(a=1.2), select=TRUE)

# how well does the model do?
print(obs_exp(dsm_boring, ~beaufort))



## 6. set up prediction data
load("data/airgrid.RData")
pred <- pred_grid[pred_grid$Stratum %in% area_of_interest, ]

# coastline for plotting
coastline <- map_data("world", c("UK", "Ireland", "Isle of Man"))

# prediction grid, three times
pp_ss <- rbind(pred, pred, pred)
pp_ss$off.set <- pp_ss$area*1000^2
pp_ss$size_class <- as.factor(c(rep(cut.labs[1], nrow(pred)),
                                rep(cut.labs[2], nrow(pred)),
                                rep(cut.labs[3], nrow(pred))))


## 7. variance estimation

# varprop
vp <- dsm_varprop(df_fs_chop, pp_ss)
# collate the abundance estimates per group and estimate variance properly
ests <- est_N_var_ss(vp, ddf_data_ss, pp_ss)

# delta method per grid cell
vg <- dsm.var.gam(df_fs_chop, split(pp_ss, 1:nrow(pp_ss)), off.set=pp_ss$off.set)
vg$pred <- unlist(vg$pred)
ests_vg <- est_N_var_ss(vg, ddf_data_ss, pp_ss)


# variance for estimated abundance model for comparison
var_gam <- dsm.var.gam(dsm_boring, list(pred), off.set = pred$area*1000^2)


## 8. actually make predictions

big_plot <- pp_ss
big_plot$Abundance <- predict(df_fs_chop, pp_ss, off.set=pp_ss$off.set)

# get abundance in each group
for(size in levels(big_plot$size_class)){
  big_plot$Abundance[big_plot$size_class==size] <-
                      big_plot$Abundance[big_plot$size_class==size] *
                      median(ddf_data_ss$save_size[ddf_data_ss$size_class==size])
}

# make combined plot predictions
comb_plot <- pred
comb_plot$Abundance_t <- 0
for(size in levels(big_plot$size_class)){
  comb_plot$Abundance_t <- comb_plot$Abundance_t +
                             big_plot$Abundance[big_plot$size_class==size] *
                      mean(ddf_data_ss$save_size[ddf_data_ss$size_class==size])
}

# make estimated abundance model predictions
boring_plot <- pred
boring_plot$Abundance <- predict(dsm_boring, pred, off.set=pred$area*1000^2)


## 9. make some plots!


# make a plot per group level
plist <- list()

for(size in levels(big_plot$size_class)){
  this_plot <- big_plot[big_plot$size_class==size,]

  plist[[size]] <- ggplot(this_plot) +
    geom_polygon(data = coastline, aes(x=long, y = lat, group = group), fill="#D2D2D2") +
    geom_tile(aes(x=Longitude, y=Latitude, fill=Abundance),
              width=0.033, height=0.066)+
    coord_quickmap(#projection="mercator",
              ylim=range(this_plot$Latitude), xlim=range(this_plot$Longitude)) +
    scale_fill_viridis()+
    geom_segment(aes(x=st.lon, y=st.lat, xend=end.lon, yend=end.lat), alpha=0.4,
                 data=segments_data_ss[segments_data_ss$size_class==size, ]) +
    geom_point(aes(x=mid.lon, y=mid.lat), alpha=0.4, size=0.25,
               data=ddf_data_ss[ddf_data_ss$size_class==size, ]) +
    ggtitle(size) +
    labs(fill="Abundance") +
    theme_minimal() + theme(axis.text=element_blank(),
                            axis.title=element_blank(),
                            plot.title=element_text(size=8),
                            legend.text=element_text(size=6),
                            legend.title=element_text(size=6),
                            legend.key.height=unit(10, "pt"),
                            legend.position="bottom")
}

# put combined group model and boring model on the same scale
# not used as we don't plot the boring model in the paper
#comb_boring_fill_scale <- c(0, max(boring_plot$Abundance, comb_plot$Abundance_t))


# make the combined plot
plist[["comb"]] <- ggplot(comb_plot) +
  geom_tile(aes(x=Longitude, y=Latitude, fill=Abundance_t),
            width=0.033, height=0.066)+
  geom_polygon(data = coastline, aes(x=long, y = lat, group = group), fill="#D2D2D2") +
  coord_quickmap(#projection="mercator",
            ylim=range(this_plot$Latitude), xlim=range(this_plot$Longitude)) +
  scale_fill_viridis()+#comb_boring_fill_scale)+
  geom_segment(aes(x=st.lon, y=st.lat, xend=end.lon, yend=end.lat), alpha=0.4,
               data=segments_data) +
  geom_point(aes(x=mid.lon, y=mid.lat), size=0.25, data=ddf_data_ss, alpha=0.4) +
  ggtitle("Combined abundance") +
  labs(fill="Total\nabundance") +
  theme_minimal() + theme(axis.text=element_blank(),
                          axis.title=element_blank(),
                          plot.title=element_text(size=8),
                          legend.text=element_text(size=6),
                          legend.title=element_text(size=6),
                          legend.key.height=unit(10, "pt"),
                          legend.position="bottom")

# test the figure that goes into the paper
#pdf("test.pdf", width=8, height=8)
#grid.arrange(grobs=list(plist[[2]],
#                        plist[[3]],
#                        plist[[1]],
#                        plist[["comb"]]), nrow=2, ncol=2)
#dev.off()

## now plot estimated abundance model
#plist[["boring"]] <- ggplot(boring_plot) +
#  geom_tile(aes(x=Longitude, y=Latitude, fill=Abundance),
#            width=0.033, height=0.066)+
#  geom_polygon(data = coastline, aes(x=long, y = lat, group = group), fill="#D2D2D2") +
#  coord_map(projection="mercator",
#            ylim=range(this_plot$Latitude), xlim=range(this_plot$Longitude)) +
#  scale_fill_viridis(comb_boring_fill_scale)+
#  geom_segment(aes(x=st.lon, y=st.lat, xend=end.lon, yend=end.lat),
#               data=segments_data_ss[segments_data_ss$size_class==size, ]) +
#  geom_point(aes(x=mid.lon, y=mid.lat),
#             data=ddf_data) +
#  ggtitle("Estimated abundance model") +
#  theme_minimal() + theme(axis.text=element_blank(),
#                          axis.title=element_blank(),
#                          legend.position="bottom")


## plot everything -- check before paper
# switching coord_map to coord_quickmap above makes the wait more
#  palatable
#grid.arrange(grobs=list(plist[[2]],
#                        plist[[3]],
#                        plist[[1]],
#                        plist[["comb"]],
#                        plist[["boring"]]), nrow=2, ncol=3)
#

## 10. save output
save(plist, df_s, df_fs_chop, vp, ests, ests_vg, dsm_boring, var_gam,
     file="hp-output.RData")

