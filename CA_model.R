library(rgrass7)
library(tools)
library(raster)
library(tidyverse)
library(tidyr)
library(sf)
library(rgeos)
library(cleangeo)
library(sp)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(readr)
library(rgdal)
library(qdapRegex)
library(mapview)
library(mapedit)
library(maptools)
library(elevatr)
library(lwgeom)
library(aod)
library(io)
library(BBmisc)
library(spatialEco)
library(caret)
library(landscapemetrics)
library(zonator)

#CA Model for TR_2 ( TR_2= [∏(NP,RP,RD,𝐒𝐌,𝐒𝐓𝐑)])
#Setting: 'Higher and lower orders of 'complete' stream channels along with neighbourhood ('B-1')
#enter the city name
city_name <- "enter CITY_NAME"
#reading the stdy area extent 
stdy_r <- raster(sprintf("xxx/stdy_r.grd",city_name)) 
plot(stdy_r)
res(stdy_r)
#reading road network of time 1(1990)
rd_netw_1_r <- raster(sprintf("xxx/rd_netw_1_r.grd",city_name)) 
plot(rd_netw_1_r)
res(rd_netw_1_r)

#reading road network of time 2 (2015)
rd_netw_2_r <- raster(sprintf("xxx/rd_netw_2_r.grd",city_name)) 
plot(rd_netw_2_r)
res(rd_netw_2_r)


#reading the stream network
strm_rd <- readOGR(dsn = sprintf("xxxx/%s",city_name), layer = "order_union")
plot(strm_rd)

#selecting only the higher order streams

strm_rd$strahler <- as.numeric(strm_rd$strahler)
strm_rd$strahler 

strm_rd$strahler_2 <- as.numeric(strm_rd$strahler_2)
strm_rd$strahler_2
#adding all strahler metrics into a single column
#ADD # BEFORE + for only drainage lines (settings D )
strm_rd$newor <- replace_na(strm_rd$strahler,0)+replace_na(strm_rd$strahler_2,0)
#only ridge (DELETE ABOVE TWO LINES AND USE THE FOLLOWING LINE FOR 'ONLY RIDGE LINES' - settings R)
#strm_rd$newor <- replace_na(strm_rd$strahler_2,0) #+ replace_na(strm_rd$strahler,0)

strm_rd$newor
range(strm_rd$newor)
strm <- st_as_sf(strm_rd)
class(strm)
max(strm$newor[])
val <- round(max(strm$newor[])/2,0)
val 
strm_1 <- strm %>% filter(newor > val)
#change 'D','R', or 'B'
strm_1 <- as(strm_1,"Spatial")
class(strm_1)
plot(strm_1)

strm_1_B <- strm_1
strm_rd_B <- strm_rd

outfile <- sprintf('xxx/%s/strm_alt_B.shp',city_name)
shapefile(strm_rd_B, outfile, overwrite = TRUE)
#only higher orders
outfile <- sprintf('xxx/%s/strm_r_B.shp',city_name)
shapefile(strm_1_B, outfile, overwrite = TRUE)
#reading rasterised streams
#all stream orders
strm_alt <- raster(sprintf("xxx/%s/rasterised/strm_alt_B.grd",city_name)) 
plot(strm_alt)
res(strm_alt)
#higher stream orders
strm_r <- raster(sprintf("xxx/%s/rasterised/strm_r_B.grd",city_name)) 
plot(strm_r)
res(strm_r)
#slope
slope <-raster(sprintf("xxx/%s/rasterised/slope.grd",city_name)) 
slope
slope <- projectRaster(slope,stdy_r)
res(slope)
res(strm_r)
res(stdy_r)
res(strm_alt)
range(slope[],na.rm = T)
slope <- tan(slope*pi/180)*100
range(slope[], na.rm = T)
replace_na(slope,Inf)
slope <- reclassify(slope, c(-1,5,5,5,20,.75,20,Inf,.1) )

fun <- function(x){x[x>=0] <- 1; return(x)}
rd_1 <- calc(rd_netw_1_r, fun)
plot(rd_1)
#write rd_1 and open in qgis to find the distance
writeRaster(rd_1, file = sprintf("E:/phd iit kgp/WORK/cellular_auto/%s/rasterised/rd_1.grd",city_name),overwrite = T)
res(rd_1)
#read road distance raster
rd_dist_1 <- raster(sprintf("E:/phd iit kgp/WORK/cellular_auto/%s/rasterised/rd_dist_1.grd",city_name)) 
plot(rd_dist_1)
res(rd_dist_1)

rd_1_co <- reclassify(rd_dist_1,
                      c(-Inf, quantile(rd_dist_1[],0.20),10,
                        quantile(rd_dist_1[],0.20),Inf,.001)) #setting the coefficients

#reading road network of time step 2
fun <- function(x){x[x>=0] <- 1; return(x)}
rd_2 <- calc(rd_netw_2_r, fun)
rd_2[]
plot(rd_2)
range(strm_r[], na.rm = T)
range(strm_alt[],na.rm = T)
(strm_alt[])
#plot(strm_r)
strm_r <- replace_na(strm_r,0)
strm_alt <- replace_na(strm_alt,0)
#important step
rd_1 <- replace_na(rd_1,0)
rd_2 <- replace_na(rd_2,0)

#---------------------
#calculate change map
ch1_2 <- rd_2-rd_1
ch1_2
#---------------------
#set weights for the neighbourhood

w <- matrix(c(0,3,0,3,9,3,0,3,0), nr = 3, nc = 3)
strm_r_alt_n <- focal(strm_alt, w)
strm_r_n <- focal(strm_r, w)
#w <- matrix(c(3,9,3), nr = 1, nc = 3)
w<-w/100
#----------------------------------
#preparation for similation
plot(rd_1)
T1 <- rd_1 #set T1 to first map
plot(T1)
print(paste("simulation starting at", Sys.time(), sep = ""))
#******************RE-ENTER START AND END YEAR****************
start_year <-1990
end_year <- 2015
total <- end_year - start_year
total
#calculate demand from rd_1 to rd_2
dfrd_1 <- as.data.frame(freq(rd_1))
dfrd_1
dfrd_2 <- as.data.frame(freq(rd_2))
dfrd_2
rddemand <- as.numeric(dfrd_1[2,2])
rddemand #
finaldemand <- as.numeric(dfrd_2[2,2])
finaldemand
andem <- (finaldemand - rddemand)/total
andem <- round(andem,0)
andem
#--------------------end calculate demand from rd_1 to rd_2
print("set demand, starting timer..")
ptm <- proc.time()
pb <- txtProgressBar(min = 0, max = total, style = 3)
writeRaster(rd_1, file = sprintf("xxx/%s/animate/%s.tif",city_name,(start_year)),overwrite = T)
#*****FOR THE TESTING TIME******

#list for stacking the rasters at each time step
s_t <- list()
for (i in 1:total){#begin running simulation
  
  urbdemand1 <- andem*1
  
  #i #
  #begin neighbourhood block
  n <- focal(T1, w=w)
  nhood <- cover(n,T1)
  nhood[] <- (1+(nhood[]-min(nhood[]))*10/(max(nhood[]-min(nhood[]))))*0.25 #re-classify so that values in range [0.25,2.75]
  plot(nhood)
  model_nhood <- nhood
  #plot(model_nhood)
  #end neighboourhood block
  #-------------------------
  #begin random block
  x <- runif(ncell(T1))
  weibull <- 1+(-log(1-(x)))*exp(1/2)
  funselect <- function(x) {x[x!=0] <- NA; return(x)}
  vacant <- calc(T1,funselect) #extract only vacant areas of T1
  random <- rd_1 # copy of the start time map
  values(random) <- weibull # assign the random values to the copy of the start time map
  model_random <- mask(random, vacant) # generating a mask to apply NAs to random raster layer
  model_random <- cover(model_random, T1) # #cover to fill the NAs back in with the original values from T1.
  model_random[] <- (1 + (model_random[]-min(model_random[]))*9/(max(model_random[])-min(model_random[])))*.025 #re-classify so that values in range [0.025,0.25]
  freq(model_random)
  #end random block
  #--------------------------------------
  #begin distance matrix for stream suitability
  #-----------------------------------------------
  T1_new <- model_nhood
  fun_T1_new <- function(x){x[x==0] <- NA; return(x)} #x>0
  T1_new <- calc(T1_new, fun_T1_new)
  
  
  strm_r_1 <- strm_r
  strm_r_1 <- strm_r_1 * (T1_new)
  freq(strm_r_1)
 
  
  
 
  #begin tansition potential calculation
 
  #replace strm_r_alt_n with strm_alt for B/D/R-3, strm_r_n for B/D/R-4, and strm_r for B/D/R-2
  model_TPP <- model_nhood*model_random*(strm_r_alt_n)*slope*(rd_1_co) 
  #replace the above line with 
  #model_TPP <- model_nhood*model_random*(rd_1_co)#for TR_1
  
 
  #coverting all 1 values of T1 into NA
  fun_T1 <- function(x){x[x==1] <- NA; return(x)}
  T1 <- calc(T1, fun_T1)
 
  model_TPP <- mask(model_TPP, T1)
  freq(model_TPP)
 
  TPP <- model_TPP
 
  max(TPP[]) #
  sort(TPP[], decreasing = T) #
  #select the highest TPP values from the TPP map
  x <- as.matrix(TPP)
  n <- urbdemand1
  x2 <- sort(-x, partial = n)
  x2h <- -sort(x2[1:n])
  ix <- which(x %in% x2h) #which positions in the matrix are the top n kept
  ix #
  y <- matrix(c(3,12,14,5,2,4,7,20,15,1), nrow = 2, ncol = 5) #
  y #
  y1 <- sort(-y, partial = 4) #
  y1 #
  y2 <- -sort(y1[1:4]) #
  y2 #
  y3 <- which(y %in% y2) #
  y3 #
  y[y3] #
  rowsT1 <- nrow(T1)
  rowsT1 #
  colsT1 <- ncol(T1)
  ro <- ix %% rowsT1
  ro #
  ro <- ifelse(ro == 0L, rowsT1, ro)
  ro #
  co <- 1 + ix %/% rowsT1
  co #
  x3 <- x[ix] #extracting the values of the top n
  x3 #
  d <- data.frame(row = ro, col = co, x = x3)
  d #
  result <- d[rev(order(d$x)),]
  result #
  #------------------------------------------
  #test for duplicates that inflate the number of cells allocated
  n #
  length(result$x) #
  difftrans <- (length(result$x)-n)
  if(difftrans>0){
    result2 <- head(result, -difftrans)#remove the duplicates from the end of the file
    result <- result2
  }
  #turn selected n values in the dataframe into a matrix and then to a raster
  x.mat <- matrix(0, rowsT1, colsT1)#create a matrix with right number of row and cols and fill with o
  x.mat #
  x.mat[cbind(result$row, result$col)] <- 1 # we want values to be 1
  max(x.mat[])
  x.mat[6,15] #
  r <- raster(x.mat)
  r[6,15] #
  extent(r) <- extent(T1)
  newdata <- r
  newdata <- mask(newdata, T1) ##generating a mask to apply NAs to final map layer
  fun_n <- function(x){x[is.na(x)] <- 1; return(x)}
  newdata <- calc(newdata, fun_n)
  
  T1 <- newdata #directly allocated all the cells at their most favourable locations acc to TPP
  
  #converting all NAs back to 1
  fun_T1 <- function(x){x[is.na(x)] <- 1; return(x)}
  T1 <- calc(T1, fun_T1)
  
  plot(T1, main = start_year+i)
  s_t[[i]] <- T1 
  #TPP
  filen <- paste("rd", (start_year+i), ".png", sep = "")
  png(filename = paste(filen),width = 1200, height = 1200, bg="white") #Plot each land use map to be able to make an animation.
 
  writeRaster(T1, file = sprintf("E:/phd iit kgp/WORK/cellular_auto/%s/animate/%s.tif",city_name,(start_year+i)),overwrite = T)
  dev.off()
  
  removeTmpFiles(h=5)
  Sys.sleep(0.1)
  setTxtProgressBar(pb,i)
  print(paste("Finished: road network simulation for ", (start_year + i), sep = ""))
  print(proc.time() - ptm)
}

close(pb)
#PLOT OF THE STARTING YEAR ROAD MAP
plot(rd_1, main = paste("initial road map - ", start_year, sep = ""))
dev.off()
rd_netw_1_r <- raster(sprintf("E:/phd iit kgp/WORK/cellular_auto/%s/rasterised/rd_netw_1_r.grd",city_name)) 
rd_1 <- replace_na(rd_1,0)

# ASSIGNING THE SIMULATED MAP TO A VARIABLE
sim99 <- T1
#STACKING THE ENDING YEAR REAL MAP AND THE SIMULATED MAP
sss <- stack(rd_2, sim99)
#ASSIGNING NAMES
name <- paste("real road network", start_year, sep ="")
name1 <- paste("real road network", end_year, sep = "")
name2 <- paste("simulated road network", end_year, sep = "")
names(sss) <- c(name1, name2)
#PLOTTING THE STACK
(plot(sss))
print(paste("simulation complete at,",Sys.time(), sep=""))


#confusion matrix and kappa index full simulated map and end year road network
sim <- as.character(sim99[])
rd <- as.character(rd_2[])
sim
conf_full <- confusionMatrix(data = factor(sim), reference = factor(rd), positive = "1")
conf_full
conf_full$overall["Kappa"]
nlevels(factor(rd))

#similarity between the real and simulated map or the AGREEMENT METRIC
#difference between the maps 
diff_rd <- rd_2-rd_1
diff_sim <- (sim99-rd_1)



plot(sim99, main = "sim")
plot(rd_2, main = "rd")
#DIFFERENCE BETWEEN THE SIMULATED MAP AND THE REAL MAP AT THE END YEAR
diff <- (sim99-rd_2)
#FALSE NEGATIVES
a <- length(sss$real.road.network2015[][sss$real.road.network2015[]==1 & sss$simulated.road.network2015[]==0])
a
#FALSE POSITIVES
b <- length(sss$real.road.network2015[][sss$real.road.network2015[]==0 & sss$simulated.road.network2015[]==1])
b
#TRUE POSITIVES
c <- length(sss$real.road.network2015[][sss$real.road.network2015[]==1 & sss$simulated.road.network2015[]==1])
c
#TRUE NEGATIVES
d <- length(sss$real.road.network2015[][sss$real.road.network2015[]==0 & sss$simulated.road.network2015[]==0])
d
#total number of cells
a+b+c+d
#number of road cells in the real road map of end year
rd_2_1s <- length(rd_2[rd_2==1])
e <- c/rd_2_1s
f <- a/rd_2_1s
g <- b/rd_2_1s
#proportion of correctly predicted road cells
e
#proportion of road cells wrongly predicted as non-road cells
f
#proportion of non-road cells wrongly predicted as road cells
g
#not similar
(10000-length(diff[diff==0]))/10000
#Compares two categorical rasters using Cohen's Kappa (d) or paired t-test statistic(s)
#sim <- as.data.frame(sim99)
#rd <- as.data.frame(rd_2)

#change maps - real and simulated
R_S <- jaccard(rd_2, sim99, x.min = 0.5, y.min = 0.5)
diff_sim <- (sim99-rd_1)
rd_1 <- raster(sprintf("xxx/%s/rasterised/real_rd_1.grd",city_name)) 
rd_2 <- raster(sprintf("xxx/%s/rasterised/real_rd_2.grd",city_name))
rd_1 <- replace_na(rd_1, 0)
rd_2 <- replace_na(rd_2, 0)
diff_rd <- rd_2-rd_1

#after rasterising data just read file using the below codes
writeRaster(diff_rd, file = sprintf("xxx/%s/rasterised/10_B_1_diff_rd.rst",city_name),overwrite = T)
writeRaster(diff_sim, file = sprintf("xxx/%s/rasterised/10_B_1_diff_sim.rst",city_name),overwrite = T)
writeRaster(rd_1, file = sprintf("xxx/%s/rasterised/10_B_1_rd_1.rst",city_name),overwrite = T)
writeRaster(rd_2, file = sprintf("xxx/%s/rasterised/10_B_1_rd_2.rst",city_name),overwrite = T)
writeRaster(sim99, file = sprintf("xxx/%s/rasterised/10_B_1_sim99.rst",city_name),overwrite = T)

plot(diff_sim, main = "sim")
plot(diff_rd, main = "rd")
#confusion matrix for the change maps - simulated and real
dif_sim <- as.character(diff_sim[])
dif_rd <- as.character(diff_rd[])
confusionMatrix(data = factor(dif_sim), reference = factor(dif_rd), positive = "1")
#the results show poor similarity between the simulated map and the map at the 'end' time. This is expected as 
#streets can intersect the territory's structuring lines at various angles rather than simply following them, 
#as simulated by the C.A. model. As a result, the resemblance of the map showing the neighbourhood of the 
#'original change map' to the 'simulated change map' is determined. 

#creating a fuzzy change map for the change map of the simulated
w2 <- matrix(c(0,0,50,0,0,0,50,50,50,0,50,50,500,50,50,0,50,50,50,0,0,0,50,0,0), nr = 5, nc = 5)
n_2 <- focal(diff_rd, w=w2)
nhood_2 <- cover(n_2,diff_rd)
fun_T5 <- function(x){x[x>0] <- 1; return(x)}
nhood_2 <- calc(nhood_2, fun_T5)
fun_T6 <- function(x){x[x<=0] <- 0; return(x)}
nhood_2 <- calc(nhood_2, fun_T6)
plot(nhood_2)
rd_mw <- as.character(nhood_2[])
#confusion matrix between the change map of the real "dif_sim" and the fuzzy change map "rd_mw" of the simulated
conf_diff_focal <- confusionMatrix(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")
conf_diff_focal
conf_diff_focal$mode
#F-meas between the real change map and fuzzy simulated change map
F_meas(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")
#sensitivity between the real change map and fuzzy simulated change map
senst_2 <- sensitivity(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")
senst_2
#stacking the fuzzy simulated change map, real change map, simulated change map
sss2 <- stack(nhood_2,diff_rd,diff_sim)
plot(sss2)
total_1 <- length(sss2$layer.2[][sss2$layer.2[]==1])  
total_1
#false positives
a1 <- length(sss2$layer.1[][sss2$layer.1[]==1 & sss2$layer.3[]==0])
a1
#false negatives
b1 <- length(sss2$layer.1[][sss2$layer.1[]==0 & sss2$layer.3[]==1])
b1
#true positives
c1 <- length(sss2$layer.1[][sss2$layer.1[]==1 & sss2$layer.3[]==1])
c1
#true negatives
d1 <- length(sss2$layer.1[][sss2$layer.1[]==0 & sss2$layer.3[]==0])
d1
sum(a1,b1,c1,d1)
#proportion of true positives by the total number of positives in fuzzy simulated change map
#which is sensitivity
c1/(c1+a1)
ref <- c1/length(sss2$layer.1[][sss2$layer.1[]==1])
ref
#specificity
specificity_2 <- d1/(d1+b1)
#balanced accuracy_BA
BA <- (specificity_2+senst_2)/2
BA
#recall
F_meas(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")

recall(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")
precision(data = factor(dif_sim), reference = factor(rd_mw), positive = "1")

rec <- c1/(c1+b1)
pre <- c1/(c1+a1)
f1 <- 2*((pre*rec)/(pre+rec))
f1

result_tab <- data.frame(name = city_name, 
                         experiment = "19_2",
                         code = "1",
                         type = "B",
                         fullnetw_accuracy = conf_full$overall["Accuracy"],
                         kappa_stat = conf_full$overall["Kappa"],
                         diff_prec = posPredValue(data = factor(dif_sim), reference = factor(rd_mw), positive = "1"),
                         diff_BA = BA,
                         diff_c = c1,
                         diff_a = a1,
                         diff_b = b1,
                         diff_d = d1,
                         R_S,
                         jacc_diff_real_sim = jaccard(diff_rd, diff_sim, x.min = 0.7, y.min = 0.7),
                         jacc_diff_neighofreal_sim = jaccard(diff_sim, nhood_2, x.min = 0.7, y.min = 0.7))
write.table(result_tab, file = 'xxx/table/result_tab_10', 
            append = T, sep = '\t', row.names = F, col.names = F)


#cleaning workspace

rm(list = ls())[grep("^x", ls())]
rm(list = ls())[grep("^*", ls())]

#repeat the code for other settings of stream channels
#For TR_1, the model_TPP <- model_nhood*model_random*(rd_1_co) replaces model_TPP <- model_nhood*model_random*(strm_r_alt_n)*slope*(rd_1_co) 
