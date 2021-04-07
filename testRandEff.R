require(tidyverse)
require(maps)
require(mapdata)
require(ggplot2)


# INLA -------------------------------------------------------
setwd('/Users/Yuki/FRA/INLAren/spde-book-files')

## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/barrier-'
)
source('R/spde-book-functions.R')

library(scales)
library(rgeos)
## High resolution maps when using map()
library(mapdata) 
## Map features, map2SpatialPolygons()
library(maptools)



# set directory and read data -----------------------------------------------------
dir1 = "/Users/Yuki/Dropbox/same"
setwd(dir1)

### for all data
temp = NULL
for(i in 1:12){
  df = read.csv(paste0("same", i, ".csv"))
  temp = rbind(temp, df)
}

summary(temp)


# make time series
times = data.frame(year = rep(min(temp$year):max(temp$year), each = length(unique(temp$month))), month = rep(1:max(temp$month)))
times = times %>% filter(month != 7) %>% filter(month != 8) %>% filter(year != 2012)
times = times %>% mutate(time = rep(1:nrow(times)))

# combine
temp = temp %>% select(year, month, lon, lat, kg) %>% filter(month != 7) %>% filter(month != 8) %>% filter(year != 2012)
unique(temp$month)
temp = left_join(temp, times, by = c("year", "month"))
summary(temp)

temp = temp %>% filter(lon < 145) %>% filter(lat < 42) %>% filter(year < 1988)
summary(temp)
catch = (temp$kg > 0) + 0 #バイナリーデータに変換
# 
# 



### for minimum data
m1 = read.csv('same1.csv')
summary(m1)
m1 = m1 %>% select(year, month, lon, lat, kg)

# 予備解析のためデータを小さくする
m1 = m1 %>% filter(year < 2012)
summary(m1)

temp = m1

# make time series
times = data.frame(year = rep(min(temp$year):max(temp$year)), month = rep(unique(temp$month)))
times = times %>% mutate(time = rep(1:nrow(times)))

temp = left_join(temp, times, by = c("year", "month"))
# temp = temp %>% filter(lon < 145) %>% filter(lat < 42)
temp = temp %>% filter(lon < 143) %>% filter(lat < 42) %>% filter(lon > 135) %>% filter(lat > 37)
summary(temp)
catch = (temp$kg > 0) + 0 #バイナリーデータに変換
summary(temp)


### for some year using all month data
m1 = NULL
for(i in 1:12){
  data = read.csv(paste0("same", i, ".csv"))
  m1 = rbind(m1, data)
}

summary(m1)
temp = m1 %>% dplyr::filter(year < 1975)
summary(temp$year)

# make time series
times = data.frame(year = rep(min(temp$year):max(temp$year)), month = rep(unique(temp$month)))
times = times %>% mutate(time = rep(1:nrow(times)))

temp = left_join(temp, times, by = c("year", "month"))
temp = temp %>% filter(lon < 145) %>% filter(lat < 42)
summary(temp)
catch = (temp$kg > 0) + 0 #バイナリーデータに変換

summary(temp)



# create mesh -----------------------------------------------------
loc = as.matrix(cbind(temp$lon, temp$lat))

mesh <- inla.mesh.2d(loc,
                     max.edge = 0.5,
                     cutoff = 0.2,
                     offset = c(1, 1))

plot(mesh)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## ----rfindex-------------------------------------------------------------
iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
                             n.group = length(unique(temp$time)))
# i = 2244

## ----apred---------------------------------------------------------------
A <- inla.spde.make.A(mesh = mesh,
                      loc = loc, group = temp$time) 
As <- inla.spde.make.A(mesh = mesh,
                      loc = loc) 


# prediction ----------------------------------------------------
# Select region 
map <- map("world", "Japan", fill = TRUE,
           col = "transparent", plot = TRUE)
IDs <- sapply(strsplit(map$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  map, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84"))
summary(map.sp)


# ## lat < 42で北海道を除去した場合
# pl.sel2 <- SpatialPolygons(list(Polygons(list(Polygon(
#   cbind(c(128, 132, 138, 144, 144, 144), # x-axis 
#         c(34,  39,  43,  43,  38,  34)), # y-axis
#   FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) #緯度経度データ

### 小さくした時
pl.sel2 <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(135, 138, 139, 139, 143, 143, 142), # x-axis
        c(37,  38.5,  40.5, 43,  43,  39,  37)), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) #緯度経度データ

summary(pl.sel2)
poly.water2 <- gDifference(pl.sel2, map.sp)
plot(pl.sel2)
plot(map.sp)
plot(poly.water2)
summary(poly.water2) #xmax = 144

tok_bor = poly.water2@polygons[[1]]@Polygons[[1]]@coords
summary(tok_bor) #xmax = 141 ここが変！=>直った？

bb_tok = poly.water2@bbox
summary(bb_tok)
x = seq(bb_tok[1, "min"] - 5, bb_tok[1, "max"] + 5, length.out = 100)
y = seq(bb_tok[2, "min"] - 5, bb_tok[2, "max"] + 5, length.out = 100)
coop = as.matrix(expand.grid(x, y))
summary(coop)
ind = point.in.polygon(coop[, 1], coop[, 2],
                       tok_bor[, 1], tok_bor[, 2])
coop = coop[which(ind == 1), ] #1516
plot(coop, asp = 1)

# check here
loc2 = NULL
for(i in 1:length(unique(temp$time))){
  loc2 = rbind(loc2, coop)
}

Ap = inla.spde.make.A(mesh = mesh, loc = loc2, group = temp$time) #error
dim(Ap) #1518, 1738
Aps = inla.spde.make.A(mesh = mesh, loc = loc2)


# stack ---------------------------------------------------------
stk = inla.stack(
  data = list(y = catch),
  A = list(A, As, 1, 1),
  effects = list(i = iset, s = spde$n.spde, intercept = rep(1, length(catch)), time = temp$time),
  tag = 'est'
)

na = rep(NA, nrow(coop)*length(unique(temp$time)))
stk.pred = inla.stack(
  data = list(y = na),
  A = list(Ap, Aps, 1, 1),
  effects = list(i = iset, s = 1:mesh$n, intercept = rep(1, nrow(coop)*length(unique(temp$time))), time = rep(1:length(unique(temp$time)), each = nrow(coop))),
  tag = 'pred'
)
joint.stk = inla.stack(stk, stk.pred)



h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))

formulae <- y ~ 0 + intercept + f(time, model = "rw1", scale.model = TRUE) + f(i, model = spde, group = i.group, control.group = list(model = 'rw1', hyper = h.spec)) + f(s, model = spde)

# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# Model fitting
res <- inla(formulae,  data = inla.stack.data(stk), 
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(stk)), 
            control.family = list(hyper = list(theta = prec.prior)), 
            control.fixed = list(expand.factor.strategy = 'inla'))

