# simple application of varprop
# minke whales from SCANS-II

library(mrds)
library(dsm)
library(ggplot2)
library(viridis)

## 1. data setup

# load minke data
load("mw_data.RData")
# from Hammond et al
trunc.ship <- 870
# truncate
ship.ddf <- ship.ddf[ship.ddf$distance <= trunc.ship, ]


## 2. fit detection functions
df_hr <- ddf(method='ds',
            dsmodel=~mcds(formula=~beaufort, key="hr"),
            data=ship.ddf, meta.data=list(width=trunc.ship))
df_hn <- ddf(method='ds',
            dsmodel=~mcds(formula=~beaufort, key="hn"),
            data=ship.ddf, meta.data=list(width=trunc.ship))
# this seems reasonable
#par(mfrow=c(1,2))
#qqplot.ddf(df_hr)
#qqplot.ddf(df_hn)
# hr definitely better by gof and aic
#> df_hn$criterion
#[1] 1025.974
#> df_hr$criterion
#[1] 1015.097
df_bf <- df_hr

## 3. fit dsm
m_ll <- dsm(count~s(x, y, k=20),# + s(depth, k=15),
            ddf.obj=df_bf, segment.data=ship.segments,
            observation.data=ship.ddf, family=tw(a=1.2))

## model with depth in didn't do too much (EDF ==1)
#summary(m_ll)
#gam.check(m_ll)
#rqgam.check(m_ll)
#source("obs_exp.R")
#obs_exp(m_ll, "beaufort")


## 4. estimate variance
vps_total <- dsm_varprop(m_ll, pred_grid)
vgs_total <- dsm.var.gam(m_ll, pred_grid, off.set=pred_grid$area*1000^2)

save(file="mw-res.RData", m_ll, pred_grid, vps_total, vgs_total, df_bf)



