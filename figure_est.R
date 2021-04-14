require(tidyverse)


setwd("/Users/Yuki/Dropbox/same")

df = NULL
m = c(1,2,3,4,5,6,9,10,11,12)
for(i in m){
  data = read.csv(paste0("est", i, ".csv")) %>% mutate(year = as.numeric(str_sub(time, 1, 4)), month = as.numeric(str_sub(time, 6,7)))
  df = rbind(df, data)
}
summary(df)

summ = df %>% na.omit() %>% group_by(year, month) %>% summarize(mean = mean(prob))

df2 = left_join(df, summ, by = c("year", "month")) %>% mutate(prob2 = exp(prob)/(1+exp(prob)), mean2 = exp(mean)/(1+exp(mean))) 
df3 = df2 %>% mutate(check = prob2-mean2) %>% filter(check > 0)
summary(df3)

t1 = rbind(df3 %>% filter(year == 1972, month > 10), df3 %>% filter(year == 1973, month < 7))
t2 = rbind(df2 %>% filter(year == 1972, month > 10) %>% na.omit() %>% filter(prob > log(0.3/0.7)), df2 %>% filter(year == 1973, month < 7) %>% na.omit() %>% filter(prob2 > 0.3))


require(ggplot2)

# with map
world_map <- map_data("world")
jap <- subset(world_map, world_map$region == "Japan")
jap_cog <- jap[jap$lat > 35 & jap$lat < 45 & jap$long > 130 & jap$long < 145, ]
pol = geom_polygon(data = jap_cog, aes(x=long, y=lat, group=group), colour="black", fill="black")
c_map = coord_map(xlim = c(134.5, 143), ylim = c(36.5, 43))

g = ggplot(df3 %>% filter(month == 1, year == 1972), aes(x = lon, y = lat, fill = prob))
g = ggplot(t2, aes(x = lon, y = lat, fill = prob2))
# r = geom_raster()
t = geom_tile()
# v = scale_fill_viridis(na.value = "transparent")
c = coord_fixed(ratio = 1)
f = facet_wrap(~ time, ncol = 8)
f = facet_wrap(~ time, ncol = 5)
labs = labs(x = "Longitude", y = "Latitude", colour = "Logit \n (encounter probability)")
labs = labs(x = "Longitude", y = "Latitude", fill = "Encounter probability")
th = theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_text(size = rel(1)),
           axis.text.y = element_text(size = rel(1)),
           axis.title.x = element_text(size = rel(1)),
           axis.title.y = element_text(size = rel(1)),
           legend.title = element_text(size = 10))
g+t+c+pol+labs+c_map+theme_bw()+scale_fill_gradientn(colours = c("blue", "cyan", "green", "yellow", "orange", "red", "darkred"))+f
