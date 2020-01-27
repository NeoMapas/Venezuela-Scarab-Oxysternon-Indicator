require(raster)

pts <- subset(Olt,!is.na(lon))
coordinates(pts) <- c("lon","lat")
proj4string(pts) <- rlim@crs

plot(rlim,legend=F)
plot(vzla,add=T)
points(pts,pch=3,cex=.3)
points(subset(Olt,pa),pch=19,col=2,cex=.8)

rpts <- distanceFromPoints(rlim,pts)
rpts <- raster(rasterize(pts,rlim),2)

plot(rpts)
points(pts,pch=3,cex=.3)

bff <- rpts<50000

plot(bff)

cbios <- rbios
values(cbios)[values(bff)<1,] <- NA

vpts <- extract(rbios,pts)
ss <- !rowSums(is.na(vpts))>0


prd1 <- raster("~/Escritorio/ME_Oxy_Normal/Oxysternon_festivum.asc")
prd2 <- raster("~/Escritorio/ME_Oxy_SWD/Oxysternon_festivum_Rasters.asc")
vls <- extract(prd1,coordinates(pts)[ss,])
vl2 <- extract(prd2,coordinates(pts)[ss,])

## Xavier Robin, Natacha Turck, Alexandre Hainard, _et al._ (2011)
##     ``pROC: an open-source package for R and S+ to analyze and compare
##     ROC curves''. _BMC Bioinformatics_, *7*, 77.  DOI:
##     10.1186/1471-2105-12-77.


require(pROC)

## CI de la sensitividad
roc1 <- roc(response=pts@data$pa[ss],predictor=vls,percent=T,ci=T,of="se",sp=seq(0,100,5))
roc2 <- roc(response=pts@data$pa[ss],predictor=vl2,percent=T,ci=T,of="se",sp=seq(0,100,5))
plot(roc1, print.auc=T, auc.polygon=F, ci.type="shape")
plot(roc2, print.auc=T, auc.polygon=F, add=T, col=2,
     print.auc.x=95, print.auc.y=95, ci.type="shape")

roc(response=pts@data$pa[ss],predictor=vls,ci=T)
roc(response=pts@data$pa[ss],predictor=vl2,ci=T)
roc.test(roc1,roc2)





df <- data.frame(species=c("background","Oxysternon festivum")[pts@data$pa[ss]+1],
           coordinates(pts)[ss,],
           vpts[ss,])           

pres <- subset(df,pts@data$pa[ss])
aus <- subset(df,!pts@data$pa[ss])

write.csv(file="~/Escritorio/presencias.csv",pres,row.names=F)
write.csv(file="~/Escritorio/ausencias.csv",aus,row.names=F)


system("mkdir ~/Escritorio/Rasters")
writeRaster(raster(rbios,1),"~/Escritorio/Rasters/arboreo.asc","ascii")
writeRaster(raster(rbios,3),"~/Escritorio/Rasters/bio01.asc","ascii")
writeRaster(raster(rbios,8),"~/Escritorio/Rasters/bio06.asc","ascii")
writeRaster(raster(rbios,18),"~/Escritorio/Rasters/bio16.asc","ascii")
writeRaster(raster(rbios,21),"~/Escritorio/Rasters/bio19.asc","ascii")

system("mkdir ~/Escritorio/Rasters_C")
writeRaster(raster(cbios,1),"~/Escritorio/Rasters_C/arboreo.asc","ascii")
writeRaster(raster(cbios,3),"~/Escritorio/Rasters_C/bio01.asc","ascii")
writeRaster(raster(cbios,8),"~/Escritorio/Rasters_C/bio06.asc","ascii")
writeRaster(raster(cbios,18),"~/Escritorio/Rasters_C/bio16.asc","ascii")
writeRaster(raster(cbios,21),"~/Escritorio/Rasters_C/bio19.asc","ascii")

system("mkdir ~/Escritorio/ME_Oxy_Normal")
system("mkdir ~/Escritorio/ME_Muestreo")
system("mkdir ~/Escritorio/ME_Oxy_SWD")
system("mkdir ~/Escritorio/ME_Oxy_SESGO")
system("mkdir ~/Escritorio/ME_Oxy_DIST")


### para Izza, variables venezuela asc
vbios <- crop(rbios,vzla)
rvzla <- rasterize(vzla,vbios)
bios.v <- (!is.na(rvzla))*vbios
values(bios.v)[is.na(values(rvzla)),] <- NA
plot(bios.v)
names(bios.v) <- names(rbios)

for (j in 1:nlayers(bios.v))
    writeRaster(raster(bios.v,j),sprintf("/media/mapoteca/IZZA/capas.asc/%s.asc",names(bios.v)[j]),"ascii")


### para otros estudiantes curso

mdr <- "/media/mapoteca/HD-E1/Mapoteca/modclim/CMIP/2_5m/hg45bi50/"
climfut <- stack(dir(mdr,full.names=T))
climfut.v <- crop(climfut,vzla)

for (j in 1:nlayers(climfut.v))
    writeRaster(raster(climfut.v,j),
                sprintf("/media/mapoteca/IZZA/futuro.asc/bio%02d.asc",j),"ascii")
