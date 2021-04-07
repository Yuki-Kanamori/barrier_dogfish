require(tidyverse)
require(ggplot2)
require(maps)


dir1 = "/Users/Yuki/Dropbox/same"

# # JS ------------------------------------------------------------
# dir_js = "/Users/Yuki/Dropbox/same/JS"
# setwd(dir_js)
# path = dir_js
# files = list.files(path)
# 
# js = NULL
# for(i in 1:length(files)){
#   data = read.csv(paste0(files[i]), fileEncoding = "CP932")
#   
#   if("Same" %in% colnames(data)){
#     data2 = data %>% select(CatchDate, Effort_tow, Same, Latitude, Longitude) %>% mutate(sp = "Same") %>% dplyr::rename(kg = Same)
#     kg = data2$kg
#     kg[is.na(kg)] = 0
#     data2$kg = kg
#   }
#   if("SameEi" %in% colnames(data)){
#     data2 = data %>% select(CatchDate, Effort_tow, SameEi, Latitude, Longitude) %>% mutate(sp = "SameEi") %>% dplyr::rename(kg = SameEi)
#     kg = data2$kg
#     kg[is.na(kg)] = 0
#     data2$kg = kg
#   }
#   js = rbind(js, data2)
# }
# unique(js$CatchDate)
# summary(js$CatchDate)
# js2 = js %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))
# summary(js2)
# unique(js2$year)
# 
# # 緯度経度が間違っている
# lat = js2 %>% filter(Latitude > 45)
# summary(lat)
# 
# 
# # setwd(dir1)
# # write.csv(js2, "js.csv", fileEncoding = "CP932")


# JS2 ------------------------------------------------------------
dir_js = "/Users/Yuki/Dropbox/same/JS2"
setwd(dir_js)
path = dir_js
files = list.files(path)

js = NULL
for(i in 1:length(files)){
  data = read.csv(paste0(files[i]), fileEncoding = "CP932")
  
  data2 = data %>% select(CatchDate, Effort_tow, Samerui, Latitude, Longitude) %>% mutate(sp = "Same") %>% dplyr::rename(kg = Samerui)
  kg = data2$kg
  kg[is.na(kg)] = 0
  data2$kg = kg
  
  js = rbind(js, data2)
}
summary(js)

js2 = js %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))
summary(js2)


# PO ------------------------------------------------------------
dir_po = "/Users/Yuki/Dropbox/same/PO"
setwd(dir_po)
path = dir_po
files = list.files(path)

po = NULL
for(i in 1:length(files)){
  data = read.csv(paste0(files[i]), fileEncoding = "CP932")
  
  if("Same" %in% colnames(data)){
    data2 = data %>% select(CatchDate, Effort_tow, Same, Latitude, Longitude) %>% mutate(sp = "Same") %>% dplyr::rename(kg = Same)
    kg = data2$kg
    kg[is.na(kg)] = 0
    data2$kg = kg
  }
  if("Aburatsunozame" %in% colnames(data)){
    data2 = data %>% select(CatchDate, Effort_tow, Aburatsunozame, Latitude, Longitude) %>% mutate(sp = "Aburatsunozame") %>% dplyr::rename(kg = Aburatsunozame)
    kg = data2$kg
    kg[is.na(kg)] = 0
    data2$kg = kg
  }
  po = rbind(po, data2)
}

unique(po$CatchDate)
summary(po$CatchDate)
po2 = po %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))
summary(po2)
unique(po2$year)

# setwd(dir1)
# write.csv(po2, "po.csv", fileEncoding = "CP932")

check = po2 %>% filter(Latitude > 44)
summary(check)

# data area -----------------------------------------------------
summary(js2)
summary(po2)

same = rbind(po2, js2) %>% dplyr::rename(lon = Longitude, lat = Latitude)
summary(same)
# 北海道も推定してみる
# same = same %>% filter(lat < 42) #うまいこと北海道のデータは除去できる

map = ggplot() + coord_fixed() + xlab("Longitude") + ylab("Latitude")
world_map = map_data("world")
region2 = subset(world_map, world_map$region == region)
local_map = map + geom_polygon(data = region2, aes(x = long, y = lat, group = group), colour = "black", fill = "white") + coord_map(xlim = c(min(same$lon)-0.2, max(same$lon)+0.2), ylim = c(min(same$lat)-0.2, max(same$lat)+0.2))
th = theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_text(size = rel(1.5)),
           axis.text.y = element_text(size = rel(1.5)),
           axis.title.x = element_text(size = rel(1.5)),
           axis.title.y = element_text(size = rel(1.5)),
           legend.title = element_text(size = 13))
p = geom_point(data = same , aes(x = lon, y = lat, colour = kg), shape = 16, size = 1)
# p = geom_point(data = same %>% filter(kg > 10000), aes(x = lon, y = lat, colour = kg), shape = 16, size = 1)
c = scale_colour_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
labs = labs(x = "Longitude", y = "Latitude", colour = "kg")
f = facet_wrap(~ year, ncol = 8)

fig = local_map+theme_bw()+th+p+c+labs+f
ggsave(filename = "map_same.pdf", plot = fig, units = "in", width = 11.69, height = 8.27)



# save data -----------------------------------------------------
for(i in 1:12){
  data = same %>% filter(month == i)
  setwd(dir1)
  file_name = paste0("same", i)
  write.csv(data, paste0(file_name, ".csv"), fileEncoding = "CP932")
}





# data check in Hokkaido ----------------------------------------
# PO ------------------------------------------------------------
dir_po = "/Users/Yuki/Dropbox/same/PO"
setwd(dir_po)
path = dir_po
files = list.files(path)

po = NULL
for(i in 1:length(files)){
  data = read.csv(paste0(files[i]), fileEncoding = "CP932")
  
  if("Same" %in% colnames(data)){
    data2 = data %>% select(CatchDate, Effort_tow, Same, Latitude, Longitude)
    kg = data2$Same
    kg[is.na(kg)] = 0
    data2$Same = kg
    
    if("Ei" %in% colnames(data)){
      ei = data$Ei
      ei[is.na(ei)] = 0
      data2 = data2 %>% mutate(Ei = ei)
    }else{
      data2$Ei = 0
    }
    
    data2$Aburatsunozame = 0
  }
  
  if("Aburatsunozame" %in% colnames(data)){
    data2 = data %>% select(CatchDate, Effort_tow, Aburatsunozame, Latitude, Longitude)
    kg = data2$Aburatsunozame
    kg[is.na(kg)] = 0
    data2$Aburatsunozame = kg
    
    data2 = data2 %>% mutate(Same = 0, Ei = 0)
  }
  po = rbind(po, data2)
}

po = po %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))

ab = po %>% filter(Aburatsunozame > 0)
summary(ab)

### check around Hoppo Ryodo
ab2 = po %>% filter(year > 1993, Longitude > 144) 
ab2 = ab2 %>% mutate(pa = ifelse(ab2$Aburatsunozame > 0, 1, 0), tag = 1)
summary(ab2)
pre = ab2 %>% dplyr::group_by(year) %>% dplyr::summarize(count = sum(pa))
rec = ab2 %>% dplyr::group_by(year) %>% dplyr::summarize(total = sum(tag))
check = left_join(pre, rec, by = "year") %>% mutate(freq = (count/total)*100)


### check around Erimo
ab3 = po %>% filter(year > 1993, Longitude < 145, Latitude > 41)
ab3 = ab3 %>% mutate(pa = ifelse(ab3$Aburatsunozame > 0, 1, 0), tag = 1)
summary(ab3)
pre2 = ab3 %>% dplyr::group_by(year) %>% dplyr::summarize(count = sum(pa))
rec2 = ab3 %>% dplyr::group_by(year) %>% dplyr::summarize(total = sum(tag))
check2 = left_join(pre2, rec2, by = "year") %>% mutate(freq = (count/total)*100)




# sukesou -----------------------------------------------------------------
require(tidyverse)
require(ggplot2)
require(maps)


dir1 = "/Users/Yuki/Dropbox/same"

# JS2 ------------------------------------------------------------
dir_js = "/Users/Yuki/Dropbox/same/JS2"
setwd(dir_js)
path = dir_js
files = list.files(path)

js = NULL
for(i in 1:length(files)){
  data = read.csv(paste0(files[i]), fileEncoding = "CP932")
  
  data2 = data %>% select(CatchDate, Effort_tow, SuketoudaraTotal, Latitude, Longitude) %>% mutate(sp = "Suke") %>% dplyr::rename(kg = SuketoudaraTotal)
  kg = data2$kg
  kg[is.na(kg)] = 0
  data2$kg = kg
  
  js = rbind(js, data2)
}
summary(js)

js2 = js %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))
summary(js2)


# PO ------------------------------------------------------------
dir_po = "/Users/Yuki/Dropbox/same/PO"
setwd(dir_po)
path = dir_po
files = list.files(path)

po = NULL
for(i in 1:length(files)){
  data = read.csv(paste0(files[i]), fileEncoding = "CP932")
  
 
  data2 = data %>% select(CatchDate, Effort_tow, Suketoudara, Latitude, Longitude) %>% mutate(sp = "Suke") %>% dplyr::rename(kg = Suketoudara)
  kg = data2$kg
  kg[is.na(kg)] = 0
  data2$kg = kg
 
  po = rbind(po, data2)
}

unique(po$CatchDate)
summary(po$CatchDate)
po2 = po %>% mutate(year = as.numeric(str_sub(CatchDate, 1, 4)), month = as.numeric(str_sub(CatchDate, 5, 6)), day = as.numeric(str_sub(CatchDate, 7, 8)))
summary(po2)
unique(po2$year)

# setwd(dir1)
# write.csv(po2, "po.csv", fileEncoding = "CP932")

check = po2 %>% filter(Latitude > 44)
summary(check)

# data area -----------------------------------------------------
summary(js2)
summary(po2)

suke = rbind(po2, js2) %>% dplyr::rename(lon = Longitude, lat = Latitude)
summary(suke)
# 北海道も推定してみる
# same = same %>% filter(lat < 42) #うまいこと北海道のデータは除去できる

map = ggplot() + coord_fixed() + xlab("Longitude") + ylab("Latitude")
world_map = map_data("world")
region2 = subset(world_map, world_map$region == region)
local_map = map + geom_polygon(data = region2, aes(x = long, y = lat, group = group), colour = "black", fill = "white") + coord_map(xlim = c(min(same$lon)-0.2, max(same$lon)+0.2), ylim = c(min(same$lat)-0.2, max(same$lat)+0.2))
th = theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_text(size = rel(1.5)),
           axis.text.y = element_text(size = rel(1.5)),
           axis.title.x = element_text(size = rel(1.5)),
           axis.title.y = element_text(size = rel(1.5)),
           legend.title = element_text(size = 13))
p = geom_point(data = same , aes(x = lon, y = lat, colour = kg), shape = 16, size = 1)
# p = geom_point(data = same %>% filter(kg > 10000), aes(x = lon, y = lat, colour = kg), shape = 16, size = 1)
c = scale_colour_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
labs = labs(x = "Longitude", y = "Latitude", colour = "kg")
f = facet_wrap(~ year, ncol = 8)

fig = local_map+theme_bw()+th+p+c+labs+f
ggsave(filename = "map_suke.pdf", plot = fig, units = "in", width = 11.69, height = 8.27)



# save data -----------------------------------------------------
for(i in 1:12){
  data = suke %>% filter(month == i)
  setwd(dir1)
  file_name = paste0("suke", i)
  write.csv(data, paste0(file_name, ".csv"), fileEncoding = "CP932")
}

