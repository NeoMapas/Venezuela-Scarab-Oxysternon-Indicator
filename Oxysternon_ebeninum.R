require(ROpenOffice)
require(gdata)
require(sp)
require(raster)

setwd("~/tmp/NeoMapas")
(load("~/Dropbox/NeoMapas/Rdata/scrbsolis.rda"))
(load("~/Dropbox/NeoMapas/Rdata/NMscrb.rda"))

cc <- grep("Oxy",colnames(scrb.solis),value=T)

ss <- subset(scrb.solis,Oxysternon.festivum>0 | Oxysternon.ebeninum>0 | Oxysternum.f..viridanum > 0)

aggregate(ss[,cc],list(ss$NM),sum)

subset(trmp.NM,vst %in% ss$vst)
subset(trmp.NM,vst %in% subset(ss,Oxysternum.f..viridanum>0)$vst)
subset(trmp.NM,vst %in% subset(ss,Oxysternon.ebeninum>0)$vst)

Oe.xy <- read.ods("~/Dropbox/NeoMapas/data/Oxysternon/LocalidadesOebeninum.ods",stringsAsFactors=F)[[2]]
paises <- shapefile("/opt/gisdata/vectorial/TMworldborders/TM_WORLD_BORDERS-0.3.shp")
COL <- shapefile("/opt/gisdata/vectorial/GADM/COL_adm1.shp")
BRA <- shapefile("/opt/gisdata/vectorial/GADM/BRA_adm1.shp")
VEN <- shapefile("/opt/gisdata/vectorial/GADM/VEN_adm1.shp")
COL0 <- shapefile("/opt/gisdata/vectorial/GADM/COL_adm0.shp")
VEN0 <- shapefile("/opt/gisdata/vectorial/GADM/VEN_adm0.shp")



coordinates(Oe.xy) <- data.frame(as.numeric(Oe.xy$long),as.numeric(Oe.xy$lat))
proj4string(Oe.xy) <- BRA@proj4string

Oe.xy@data$ADM <- over(Oe.xy,VEN)$NAME_1
Oe.xy@data$ADM <- ifelse(is.na(Oe.xy@data$ADM),over(Oe.xy,COL)$NAME_1,paste(Oe.xy@data$ADM,"(VEN)"))
Oe.xy@data$ADM <- ifelse(is.na(Oe.xy@data$ADM),paste(over(Oe.xy,BRA)$NAME_1,"(BRA)"),Oe.xy@data$ADM)
vs <- rgb(c(184,174)/255,
          c(213,232)/255,
          c(138,84)/255,
          .25)
as <- rgb(c(231,224)/255,
          c(232,225)/255,
          c(84,150)/255,
          .25)
bs <- rgb(c(150,107)/255,
          c(206,184)/255,
          c(225,167)/255,
          .25)


e0 <- extent(Oe.xy)

png(file="Oebeninum_Fig1.png")
layout(1)
plot(Oe.xy,pch=NA)
plot(VEN,col=bs,add=T,border="grey77")
plot(COL,col=as,add=T,border="grey77")
plot(BRA,col=vs,add=T,border="grey77")
abline(h=seq(-2,6,by=2),lty=3)
abline(v=seq(-70,-60,by=2),lty=3)
plot(VEN0,add=T,lwd=3,border="grey58")
plot(COL0,add=T,lwd=3,border="grey58")
points(Oe.xy,cex=1.5,col="slateblue4",
       pch=c(1,19)[(Oe.xy@data$source %in% "NM")+1],lwd=3)

text(c(-65.93336, -63.40620, -62.59501, -65.00000, -68.64053, -68.58533),
     c(3.3622612,  5.3382325,  2.8488198, -0.4007108,  2.8220844,  5.3537914),
     c("Amazonas\n(VEN)","Bolívar\n(VEN)", "Roraima\n(BRA)","Amazonas\n(BRA)","Guainía\n(COL)","Vichada\n(COL)"),
     col="grey38",font=2)
axis(1,seq(-70,-60,by=2),sprintf("%s° W",abs(seq(-70,-60,by=2))))
axis(2,seq(-2,6,by=2),sprintf("%s° %s",abs(seq(-2,6,by=2)),c("S","","N","N","N")))
box()

## aquí me falta un registro del Brasil
bios <- stack(dir("~/SIG/mapas/Venezuela/WC_2.5min/","bio[0-9]+.tif",full.names=T))
qry <- extract(bios,Oe.xy)
m <- apply(qry,2,mean,na.rm=T)
sg <- apply(qry,2,sd,na.rm=T)

stbio <- t((t(qry)-m)/sg)

idx <- cellFromXY(bios,Oe.xy)

dts <- cbind(Oe.xy@data[!duplicated(idx),],
             stbio[!duplicated(idx),])

      
      

a1 <- subset(dts,source != "NM" & !is.na(bio01))[,-(1:4)]
a2 <- subset(dts,!is.na(bio01))[,-(1:4)]
require(Hmisc)
v <- varclus(as.matrix(a1), similarity="spear")

##computacionalmente singular
##v <- varclus(as.matrix(a1), similarity="pearson")
##v <- varclus(as.matrix(a1), similarity="hoeffding")
plot(v)

slc <- !duplicated(cutree(v$hclust,h=0.3))
colnames(a1)[slc]

mahalanobis(a1[,slc], center=apply(a1[,slc],2,mean),cov=cov(a1[,slc]))
D2 <- mahalanobis(a2[,slc], center=apply(a1[,slc],2,mean),cov=cov(a1[,slc]))

evl <- 1 - pchisq(D2, sum(slc) - 1)

h0 <- hclust(dist(a2[,slc]))

png(file="Oebeninum_Fig2.png")
layout(1:2)
par(mar=c(0,5,0,3))
plot(h0,label=subset(dts,!is.na(bio01))$ADM,xlab="",sub="",main="",hang=-2)
par(mar=c(1,5,1,3))
plot(1*D2[h0$order],pch=c(1,19)[(subset(dts,!is.na(bio01))$source %in% "NM")[h0$order] + 1],ylab=expression(D^2),xlab="",axes=F,ylim=c(90,0))
axis(2,seq(0,80,by=20),seq(0,80,by=20))
axis(3,1:10,rep("",10))
##axis(1,1:10,subset(dts,!is.na(bio01))$ADM,las=2)
box()
abline(h=10,lty=2,col=2)
dev.off()

round(evl,3)

layout(matrix(1:6,ncol=2))
for (i in 1:6)
 plot(qry[!duplicated(idx),colnames(a1)[slc]][,i])


##text(VEN,"NAME_1")
##text(COL,"NAME_1") ## Vichada, Guainía
##text(BRA,"NAME_1") ## Roraima, Amazonas

plot(paises,add=T)

vs <- rgb(c(184,174)/255,
          c(213,232)/255,
          c(138,84)/255)




