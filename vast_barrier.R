# Packages
require(TMB)
require(VAST)
require(tidyverse)
require(maps)
require(mapdata)
require(ggplot2)

# set directory and read data -----------------------------------------------------
dir1 = "/Users/Yuki/Dropbox/same"
setwd(dir1)

month = c(1,2,3,4,5,6,9,10,11,12)

m = 1
m1 = read.csv(paste0("same", m, ".csv"))
summary(m1)
m1 = m1 %>% select(year, month, lon, lat, kg, Effort_tow)

# 予備解析のためデータを小さくする
m1 = m1 %>% filter(year < 1982)
summary(m1)

temp = m1

# make time series
times = data.frame(year = rep(min(temp$year):max(temp$year)), month = rep(unique(temp$month)))
times = times %>% mutate(time = rep(1:nrow(times)))

temp = left_join(temp, times, by = c("year", "month"))
# temp = temp %>% filter(lon < 145) %>% filter(lat < 42)
temp = temp %>% filter(lon < 143) %>% filter(lat < 42) %>% filter(lon > 135) %>% filter(lat > 37)
summary(temp)
# catch = (temp$kg > 0) + 0 #バイナリーデータに変換
# summary(temp)

df = temp %>% mutate(cpue = kg/Effort_tow)

# 1. Settings ------------------------------------------------------
# 1.1 Version for cpp code
Version = get_latest_version(package = "VAST")

# 1.2 Spatial settings
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
Kmeans_Config = list("randomseed" = 1, "nstart" = 100, "iter.max" = 1000)
grid_size_km = 25
temp = temp %>% mutate(tag = paste(lon, lat, sep = "_")) 
n_x = length(unique(temp$tag))

# 1.3 Model settings
FieldConfig = c(Omega1 = 1, Epsilon1 = 1, Omega2 = 1, Epsilon2 = 1) #factor analysis
RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0) #0: fixed, 1: independent, 2:RW, 3:constant, 4:AR
OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0) #overdispersion
ObsModel = c(PosDist = 2, Link = 0)
Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0, 
            Calculate_effective_area = 1, Calculate_Cov_SE = 0, 
            Calculate_Synchrony = 0, Calculate_Coherence = 0)

# 1.4 Stratification for results
strata.limits = data.frame('STRATA'="All_areas")

# 1.5 Derived objects
Region = "others"

# 1.6 Save settings
dirname = dir1
DateFile = paste0(dirname,'/same', m, Sys.Date())
dir.create(DateFile)
Record = list(Version = Version, Method = Method, grid_size_km = grid_size_km, n_x = n_x, 
              FieldConfig = FieldConfig, RhoConfig = RhoConfig, OverdispersionConfig = OverdispersionConfig, 
              ObsModel = ObsModel, Kmeans_Config = Kmeans_Config, Region = Region,
              strata.limits = strata.limits) 
setwd(dir = DateFile)
save(Record, file = file.path(DateFile, "Record.RData")) 
capture.output(Record, file = paste0(DateFile, "/Record.txt"))


# 2. Prepare the data ----------------------------------------------
# 2.1 Data-frame
head(df)
Data_Geostat = df %>% select(year, lon, lat, cpue) %>% rename(Year = year, Lon = lon, Lat = lat, Catch_KG = cpue)

# 2.2 Extrapolation grid
Extrapolation_List = FishStatsUtils::make_extrapolation_info(
  Regio = Region, #zone range in Japan is 51:56
  strata.limits = strata.limits, 
  observations_LL = Data_Geostat[, c("Lat", "Lon")]
)

# 2.3 derived objects for spatio-temporal estimation
knot_method = "samples"
LON_intensity = Lon_i = Data_Geostat[, "Lon"]
LAT_intensity = Lat_i = Data_Geostat[, "Lat"]

loc_i = project_coordinates( X=Lon_i, Y=Lat_i, projargs=Extrapolation_List$projargs ) #UTM変換
loc_intensity = project_coordinates( X=LON_intensity, Y=LAT_intensity, projargs=Extrapolation_List$projargs ) #上と何が違う？
colnames(loc_i) = colnames(loc_intensity) = c("E_km", "N_km")
Grid_bounds = grid_size_km * apply(Extrapolation_List$Data_Extrap[,c('E_km','N_km')]/grid_size_km, MARGIN=2, FUN=function(vec){trunc(range(vec))+c(0,1)}) #範囲の最大と最小を出している？2行2列のデータが返される

# Calculate k-means centroids
make_kmeans <-
  function( n_x, loc_orig, nstart=100, randomseed=1, iter.max=1000, DirPath=paste0(getwd(),"/"),
            Save_Results=TRUE, kmeans_purpose="spatial", backwards_compatible_kmeans=FALSE ){
    
    # get old seed
    oldseed = ceiling(runif(1,min=1,max=1e6))
    # fix new seed
    if( !is.null(randomseed) ) set.seed( round(randomseed) )
    old.options <- options()
    options( "warn" = -1 )
    on.exit( options(old.options) )
    
    # warnings
    if( paste0("Kmeans-",n_x,".RData") %in% list.files(DirPath) ){
      warning("Found `Kmeans-",n_x,".RData` in directory; this might indicate a desire to import previously constructed k-means under deprecated naming conventions")
    }
    
    if(kmeans_purpose=="spatial"){
      tmpfile <- paste0("Kmeans_knots-",n_x,".RData")
    } else if(kmeans_purpose=="extrapolation"){
      ## n_x is really max_cells in this case
      tmpfile <- paste0("Kmeans_extrapolation-",n_x,".RData")
    } else {
      stop("Invalid kmeans_purpose for make_kmeans:", kmeans_purpose)
    }
    
    # Backwards compatibility
    if( backwards_compatible_kmeans==TRUE ){
      if( identical(formalArgs(RNGkind), c("kind","normal.kind","sample.kind")) ){
        RNGkind_orig = RNGkind()
        on.exit( RNGkind(kind=RNGkind_orig[1], normal.kind=RNGkind_orig[2], sample.kind=RNGkind_orig[3]), add=TRUE )
        RNGkind( sample.kind="Rounding" )
      }else if( !identical(formalArgs(RNGkind), c("kind","normal.kind")) ){
        stop("Assumptions about `RNGkind` are not met within `make_kmeans`; please report problem to package developers")
      }
    }
    
    # Calculate knots for SPDE mesh
    if( length(unique(paste(loc_orig[,1],loc_orig[,2],sep="_")))<=n_x ){
      # If number of knots is less than number of sample locations
      Kmeans = NULL
      Kmeans[["centers"]] = unique( loc_orig )
      Kmeans[["cluster"]] = RANN::nn2( data=Kmeans[["centers"]], query=loc_orig, k=1)$nn.idx[,1]
      message( "n_x greater than or equal to n_unique so no calculation necessary" )
    }else{
      if(tmpfile  %in% list.files(DirPath) ){
        # If previously saved knots are available
        load( file=paste0(DirPath,"/",tmpfile))
        message( "Loaded from ", DirPath, "/", tmpfile)
      }else{
        # Multiple runs to find optimal knots
        message("Using ", nstart, " iterations to find optimal ",
                ifelse(kmeans_purpose=="extrapolation", "extrapolation grid", "spatial knot"),
                " placement because no saved file found...")
        Kmeans = list( "tot.withinss"=Inf )
        tries <- 1
        for(i in 1:nstart){
          Tmp = stats::kmeans( x=loc_orig, centers=n_x, iter.max=iter.max, nstart=1, trace=0)
          if(i==1) message("Iter=1: Current=", round(Tmp$tot.withinss,0))
          if(i!=1 & i!=nstart & i %% 10 ==0)
            message( 'Iter=',i,': Current=',round(Kmeans$tot.withinss,0),' Proposed=',round(Tmp$tot.withinss,0) )#,' Time=',round(Time,4)) )
          if( Tmp$tot.withinss < Kmeans$tot.withinss ){
            Kmeans = Tmp
            tries <- i # which iteration was the optimal
          }
        }
        message("Iter=", nstart, ': Final=', round(Kmeans$tot.withinss,0), " after ", tries, " iterations")
        if(Save_Results==TRUE){
          save( Kmeans, file=paste0(DirPath, tmpfile))
          message( "Results saved to ", DirPath, tmpfile, "\n for subsequent runs by default (delete it to override)")
        }
      }
    }
    
    # fix to old seed
    if( !is.null(randomseed) ) set.seed( oldseed )
    
    # Return stuff
    Return = list("centers"=Kmeans[["centers"]], "cluster"=Kmeans[["cluster"]] )
    return( Return )
  }

Kmeans = make_kmeans(n_x=n_x, loc_orig=loc_intensity[,c("E_km", "N_km")], randomseed=1, nstart=100, DirPath=paste0(getwd(),"/"), Save_Results=TRUE, backwards_compatible_kmeans=FALSE)
NN_i = RANN::nn2( data=Kmeans[["centers"]], query=loc_i, k=1)$nn.idx[,1] #RANNはk近傍法関係のパッケージでリスト（nn.idxとnn.dists）が返される．行の数はDGの行数と同じ．それぞれの地点がどのknotになるのかを作成している？
length(unique(NN_i))

# Calculate grid for 2D AR1 process
loc_grid = expand.grid( 'E_km'=seq(Grid_bounds[1,1],Grid_bounds[2,1],by=grid_size_km), 'N_km'=seq(Grid_bounds[1,2],Grid_bounds[2,2],by=grid_size_km) ) #300行2列のdataframe（E_km, N_km）
Which = sort(unique(RANN::nn2(data=loc_grid, query=Extrapolation_List$Data_Extrap[which(Extrapolation_List$Area_km2_x>0),c('E_km','N_km')], k=1)$nn.idx[,1])) #loc_gridに対してknotを作成し，ユニークなものを抽出？63行
loc_grid = loc_grid[Which,] #whichの行だけ抽出
grid_num = RANN::nn2( data=loc_grid, query=loc_i, k=1)$nn.idx[,1] #もう一回knotを作成？これはよくわからん．行の数はDFの行数と同じ

# Calc design matrix and areas
knot_i = NN_i
loc_x = Kmeans[["centers"]]

# fine_scale==FALSEの時
# loc_g = loc_x

# fine_scale==TRUEの時
loc_g = Extrapolation_List$Data_Extrap[ which(Extrapolation_List$Area_km2_x>0), c('E_km','N_km') ]

# Convert loc_x back to location in lat-long coordinates latlon_x

origargs = "+proj=longlat +datum=WGS84"
latlon_x = project_coordinates( X=loc_x[,"E_km"], Y=loc_x[,"N_km"], projargs=origargs, origargs=Extrapolation_List$projargs )[,c("Y","X")] #UTMを緯度経度に戻す
colnames(latlon_x) = c("Lat", "Lon")

# Convert loc_g back to location in lat-long coordinates latlon_g
latlon_g = project_coordinates( X=loc_g[,"E_km"], Y=loc_g[,"N_km"], projargs=origargs, origargs=Extrapolation_List$projargs )[,c("Y","X")]
colnames(latlon_g) = c("Lat", "Lon")

# Bundle lat-lon
latlon_i = cbind( 'Lat'=Lat_i, 'Lon'=Lon_i )

### make_mesh
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

# North Japan -----------------------------------------------------------
# Select region 
# map <- map("world", "Japan", fill = TRUE,
#            col = "transparent", plot = TRUE)Russia
map <- map("world", c("Japan", "Russia"), fill = TRUE,
           col = "transparent", plot = TRUE)
IDs <- sapply(strsplit(map$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  map, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84")) #緯度経度データ
summary(map.sp)



# make a polygon ------------------------------------------------------------

### 小さくした時
pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(135, 138, 139, 139, 143, 143, 142), # x-axis
        c(37,  38.5,  40.5, 43,  43,  39,  37)), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) #緯度経度データ

summary(pl.sel)
poly.water <- gDifference(pl.sel, map.sp) #緯度経度
summary(poly.water)
plot(pl.sel)
plot(map.sp)
plot(poly.water)

## ------------------------------------------------------------------------
mesh.not <- inla.mesh.2d(boundary = poly.water, max.edge = 0.5,
                         cutoff = 0.2)
plot(mesh.not)
summary(mesh.not) #緯度経度


## ----label = "plot-barr-mesh1", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 10, heigh = 4.5, width = '97%', fig.cap = "The left plot shows the polygon for land in grey and the manually constructed polygon for our study area in light blue. The right plot shows the simple mesh, constructed only in the water."----
par(mfrow = c(1, 1), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.7, 0), las = 1)
par(mar = c(0, 0, 0, 0))

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(map.sp, add = TRUE, col = alpha(gray(0.9), 0.5))

plot(pl.sel, asp = 1)
plot(map.sp, add = TRUE, col = gray(0.9))
plot(mesh.not, add = TRUE) #ここが変！=>修正済

## ------------------------------------------------------------------------
mesh <- inla.mesh.2d(boundary = poly.water,
                     max.edge = 0.5,
                     cutoff = 0.2,
                     offset = c(1, 1))
plot(mesh)

anisotropic_mesh = mesh

## ------------------------------------------------------------------------
water.tri = inla.over_sp_mesh(poly.water, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
# ignore.CRS: whether to ignore the coordinate system information in x and y (default FALSE)
summary(water.tri)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, 
                                    barrier.triangles = barrier.tri)
plot(poly.barrier)

## ----label = "plot-barr-mesh2", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 6, heigh = 4.5, width = '97%', fig.cap = "The mesh constructed both over water and land. The grey region is the original land map. The inner red outline marks the coastline barrier."----

plot(mesh, lwd = 0.5, add = FALSE)
plot(pl.sel, add = TRUE)
plot(map.sp, add = TRUE, col = gray(.9))
plot(poly.barrier, border = "red", add = TRUE)


# ---------------------------------------------------------------
# spde ----------------------------------------------------------
# ---------------------------------------------------------------
range <- 200
barrier.model <- inla.barrier.pcmatern(mesh, 
                                       barrier.triangles = barrier.tri)
anisotropic_spde = barrier.model


# Pre-processing in R for anisotropy
Dset = 1:2
# Triangle info
TV = anisotropic_mesh$graph$tv       # Triangle to vertex indexing
V0 = anisotropic_mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = anisotropic_mesh$loc[TV[,2],Dset]
V2 = anisotropic_mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0

barrier_finite_elements = INLA:::inla.barrier.fem(mesh=mesh,
                                                  barrier.triangles=barrier.tri)
barrier_list = list(C0 = barrier_finite_elements$C[[1]],
                    C1 = barrier_finite_elements$C[[2]],
                    D0 = barrier_finite_elements$D[[1]],
                    D1 = barrier_finite_elements$D[[2]],
                    I = barrier_finite_elements$I )

# Calculate Areas #
crossprod_fn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = crossprod_fn( E0[i,],E1[i,] )/2   # T = area of each triangle

################
# Add the isotropic SPDE mesh for spherical or 2D projection, depending upon `Method` input
################

# Mesh and SPDE for different inputs
loc_isotropic_mesh = loc_x
isotropic_mesh = anisotropic_mesh
# isotropic_spde = INLA::inla.spde2.matern(isotropic_mesh, alpha=2)
isotropic_spde = anisotropic_spde

MeshList = list("loc_x"=loc_x, "loc_isotropic_mesh"=loc_isotropic_mesh, "isotropic_mesh"=isotropic_mesh,
                "isotropic_spde"=isotropic_spde, "anisotropic_mesh"=anisotropic_mesh, "anisotropic_spde"=anisotropic_spde,
                "Tri_Area"=Tri_Area, "TV"=TV, "E0"=E0, "E1"=E1, "E2"=E2,
                #"anisotropic_mesh_triangles_over_land"=anisotropic_mesh_triangles_over_land, 
                "barrier_list"=barrier_list )

n_s = switch( tolower(Method), "mesh"=MeshList$anisotropic_spde$f$n, "grid"=nrow(loc_x),
              "spherical_mesh"=MeshList$isotropic_spde$n.spde, "stream_network"=nrow(loc_x), "barrier"=MeshList$anisotropic_spde$n.spde,  ) #toloweは大文字を小文字に変換（逆はtoupper）．n.spdeはbarrierにした場合はbarrier.model$f$nかな？

# Make matrices for 2D AR1 process
Dist_grid = dist(loc_grid, diag=TRUE, upper=TRUE)
M0 = as( ifelse(as.matrix(Dist_grid)==0, 1, 0), "dgTMatrix" )
M1 = as( ifelse(as.matrix(Dist_grid)==grid_size_km, 1, 0), "dgTMatrix" )
M2 = as( ifelse(as.matrix(Dist_grid)==sqrt(2)*grid_size_km, 1, 0), "dgTMatrix" )
if( Method=="Spherical_mesh" ) GridList = list("M0"=M0, "M1"=M1, "M2"=M2, "grid_size_km"=grid_size_LL)
if( Method %in% c("Mesh","Grid","Stream_network","Barrier") ) GridList = list("M0"=M0, "M1"=M1, "M2"=M2, "grid_size_km"=grid_size_km)

# Make projection matrices
# fine_scale==FALSEの時
# A_is = matrix(0, nrow=nrow(loc_i), ncol=116)
# A_is[ cbind(1:nrow(loc_i),knot_i) ] = 1
# A_is = as( A_is, "dgTMatrix" )
# A_gs = as( diag(n_x), "dgTMatrix" )
# 
# fine_scale==TRUEの時
A_is = INLA::inla.spde.make.A( MeshList$anisotropic_mesh, loc=as.matrix(loc_i) )
if( class(A_is)=="dgCMatrix" ) A_is = as( A_is, "dgTMatrix" )
A_gs = INLA::inla.spde.make.A( MeshList$anisotropic_mesh, loc=as.matrix(loc_g) )
if( class(A_gs)=="dgCMatrix" ) A_gs = as( A_gs, "dgTMatrix" )
Check_i = apply( A_is, MARGIN=1, FUN=function(vec){sum(vec>0)})
Check_g = apply( A_is, MARGIN=1, FUN=function(vec){sum(vec>0)})

# Calculate areas
if( Method != "Stream_network" ){
  # ここ
  PolygonList = Calc_Polygon_Areas_and_Polygons_Fn( loc_x=loc_x, Data_Extrap=Extrapolation_List[["Data_Extrap"]], a_el=Extrapolation_List[["a_el"]])
  if( fine_scale==FALSE ){
    # a_gl = PolygonList[["a_xl"]] 本当はこっち
    a_gl = as.matrix(Extrapolation_List[["a_el"]][ which(Extrapolation_List$Area_km2_x>0), ])
  }
  if( fine_scale==TRUE ){
    a_gl = as.matrix(Extrapolation_List[["a_el"]][ which(Extrapolation_List$Area_km2_x>0), ])
  }
}else{
  PolygonList = NULL
  dist_inp = Network_sz_LL[,"dist_s"]
  dist_inp[which(is.infinite(dist_inp))] <- 0
  a_gl = matrix(dist_inp, nrow=n_x)
}

# Moving これ何？
if( fine_scale==TRUE | Method=="Stream_network" ){
  g_e = rep(NA, length(Extrapolation_List[["Area_km2_x"]]))
  g_e[ which(Extrapolation_List[["Area_km2_x"]]>0) ] = 1:length(which(Extrapolation_List[["Area_km2_x"]]>0))
}else{
  # ここ
  g_e = PolygonList$NN_Extrap$nn.idx[,1]
  g_e[ which(Extrapolation_List[["Area_km2_x"]]==0) ] = NA
}

fine_scale = TRUE
Spatial_List = list( "fine_scale"=fine_scale, "A_is"=A_is, "A_gs"=A_gs, "n_x"=n_x, "n_s"=n_s, "n_g"=nrow(a_gl), "n_i"=nrow(loc_i),
                     "MeshList"=MeshList, "GridList"=GridList, "a_gl"=a_gl, "a_xl"=a_gl, "Kmeans"=Kmeans, "knot_i"=knot_i,
                     "loc_i"=as.matrix(loc_i), "loc_x"=as.matrix(loc_x), "loc_g"=as.matrix(loc_g), "g_e"=g_e,
                     "Method"=Method, "PolygonList"=PolygonList, "NN_Extrap"=PolygonList$NN_Extrap, "knot_method"=knot_method,
                     "latlon_x"=latlon_x, "latlon_g"=latlon_g, "latlon_i"=latlon_i )

# 2.4 save Data_Geostat
Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List[["knot_i"]], zone = Extrapolation_List[["zone"]])
write.csv(Data_Geostat, "Data_Geostat.csv")

# 3. Build and run model ----------------------------------------
# L809でエラー発生
TmbData = make_data(
  Version = Version,
  FieldConfig = FieldConfig, 
  OverdispersionConfig = OverdispersionConfig, 
  RhoConfig = RhoConfig, 
  ObsModel = ObsModel, 
  c_iz = rep(0, nrow(Data_Geostat)), 
  b_i = Data_Geostat[, 'Catch_KG'], 
  a_i = rep(1, nrow(Data_Geostat)), #cpue
  #a_i = Data_Geostat[, 'AreaSwept_km2'], # catch and effort
  s_i = Data_Geostat[, 'knot_i'] - 1,
  t_i = Data_Geostat[, 'Year'], 
  spatial_list = Spatial_List, 
  Options = Options,
  Aniso = TRUE
)

TmbList = VAST::make_model(TmbData = TmbData,
                           RunDir = DateFile,
                           Version = Version,
                           RhoConfig = RhoConfig,
                           loc_x = Spatial_List$loc_x,
                           Method = Spatial_List$Method)

Obj = TmbList[["Obj"]]
Opt = TMBhelper::fit_tmb(obj = Obj, 
                         lower = TmbList[["Lower"]], 
                         upper = TmbList[["Upper"]],
                         getsd = TRUE, 
                         savedir = DateFile, 
                         bias.correct = TRUE)
(check = VAST::check_fit(Opt))

Report = Obj$report()
Save = list("Opt" = Opt, 
            "Report" = Report, 
            "ParHat" = Obj$env$parList(Opt$par),
            "TmbData" = TmbData)
save(Save, file = paste0(DateFile,"/Save.RData"))


# 4. Figures -------------------------------------------------------
# 4.1 Plot data
plot_data(Extrapolation_List = Extrapolation_List, Spatial_List = Spatial_List, Data_Geostat = Data_Geostat, PlotDir = DateFile)

# 4.2 Convergence
pander::pandoc.table(Opt$diagnostics[, c('Param','Lower','MLE','Upper','final_gradient')])

# 4.3 Diagnostics for encounter-probability component
Enc_prob = plot_encounter_diagnostic(Report = Report, Data_Geostat = Data_Geostat, DirName = DateFile)

# 4.4 Diagnostics for positive-catch-rate component
Q = plot_quantile_diagnostic(TmbData = TmbData, 
                             Report = Report, 
                             FileName_PP = "Posterior_Predictive",
                             FileName_Phist = "Posterior_Predictive-Histogram", 
                             FileName_QQ = "Q-Q_plot", 
                             FileName_Qhist = "Q-Q_hist", 
                             DateFile = DateFile )

# 4.5 Diagnostics for plotting residuals on a map
MapDetails_List = make_map_info("Region" = Region, 
                                "spatial_list" = Spatial_List, 
                                "Extrapolation_List" = Extrapolation_List)
Year_Set = seq(min(Data_Geostat[,'Year']), max(Data_Geostat[,'Year']))
Years2Include = which(Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

plot_residuals(Lat_i = Data_Geostat[,'Lat'], 
               Lon_i = Data_Geostat[,'Lon'], 
               TmbData = TmbData, 
               Report = Report, 
               Q = Q, 
               savedir = DateFile, 
               spatial_list = Spatial_List, #
               extrapolation_list = Extrapolation_List, #
               MappingDetails = MapDetails_List[["MappingDetails"]], 
               PlotDF = MapDetails_List[["PlotDF"]], 
               MapSizeRatio = MapDetails_List[["MapSizeRatio"]], 
               Xlim = MapDetails_List[["Xlim"]], 
               Ylim = MapDetails_List[["Ylim"]], 
               FileName = DateFile, 
               Year_Set = Year_Set, 
               Years2Include = Years2Include, 
               Rotate = MapDetails_List[["Rotate"]], 
               Cex = MapDetails_List[["Cex"]], 
               Legend = MapDetails_List[["Legend"]], 
               zone = MapDetails_List[["Zone"]], 
               mar = c(0,0,2,0), 
               oma = c(3.5,3.5,0,0), 
               cex = 1.8)

# 4.6 Direction of "geometric anisotropy"
plot_anisotropy(FileName = paste0(DateFile,"Aniso.png"), 
                Report = Report, 
                TmbData = TmbData)

# 4.7 Density surface for each year
Dens_xt = plot_maps(plot_set = c(3), 
                    MappingDetails = MapDetails_List[["MappingDetails"]], 
                    Report = Report, 
                    Sdreport = Opt$SD, 
                    PlotDF = MapDetails_List[["PlotDF"]], 
                    MapSizeRatio = MapDetails_List[["MapSizeRatio"]], 
                    Xlim = MapDetails_List[["Xlim"]], 
                    Ylim = MapDetails_List[["Ylim"]], 
                    FileName = DateFile, 
                    Year_Set = Year_Set, 
                    Years2Include = Years2Include, 
                    Rotate = MapDetails_List[["Rotate"]],
                    Cex = MapDetails_List[["Cex"]], 
                    Legend = MapDetails_List[["Legend"]], 
                    zone = MapDetails_List[["Zone"]], 
                    mar = c(0,0,2,0), 
                    oma = c(3.5,3.5,0,0), 
                    cex = 1.8, 
                    plot_legend_fig = FALSE)

Dens_DF = cbind("Density" = as.vector(Dens_xt), 
                "Year" = Year_Set[col(Dens_xt)], 
                "E_km" = Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], 
                "N_km" = Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'])

pander::pandoc.table(Dens_DF[1:6,], digits=3)

# 4.8 Index of abundance
Index = plot_biomass_index(DirName = DateFile, 
                           TmbData = TmbData, 
                           Sdreport = Opt[["SD"]], 
                           Year_Set = Year_Set, 
                           Years2Include = Years2Include, 
                           use_biascorr=TRUE )
pander::pandoc.table(Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] ) 

# 4.9 Center of gravity and range expansion/contraction
plot_range_index(Report = Report, 
                 TmbData = TmbData, 
                 Sdreport = Opt[["SD"]], 
                 Znames = colnames(TmbData$Z_xm), 
                 PlotDir = DateFile, 
                 Year_Set = Year_Set)

