### R code from vignette source '~/Dropbox/NeoMapas/doc/902_AnalisisMorfometriaEscarabajos/Documento2_Oxysternonfestivum.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: SIG
###################################################

if (!exists("rbios")) {
  load("~/NeoMapas/Rdata/SIG.rda")
  rbios <- stack(c("~/NeoMapas/data/gisdata/Ofestivum/alt.tif",sprintf("~/NeoMapas/data/gisdata/Ofestivum/bio%s.tif",1:19)))
  values(rbios)[values(rbios)>55000] <- NA
  rarbo <- raster(sprintf("~/NeoMapas/data/gisdata/Ofestivum/arboreo.tif"))
  values(rarbo)[values(rarbo)==0] <- NA
  values(rarbo)[values(rarbo)==254] <- 0
  values(rarbo)[values(rarbo)==255] <- 5
  r1 <- aggregate(rarbo,5,fun=mean,na.rm=TRUE)
  arboreo <- raster::resample(r1, raster(rbios,1), method="bilinear")
  rbios <- stack(arboreo,rbios)
  Olit <- read.ods(file="~/NeoMapas/data/Oxysternon/LocalidadesOfestivum.ods")
  Oshp <- shapefile(file="~/NeoMapas/data/Oxysternon/Ofestivum.shp")
  ##bbox -73.3755402 -46.30876 -3.470826 13
  rlim <- rasterize(Oshp,arboreo)
  sbios <- rbios*rlim
  sbios@layernames <- c("arboreo","alt",sprintf("bio%02d",1:19))
  rbios@layernames <- c("arboreo","alt",sprintf("bio%02d",1:19))
}



###################################################
### code chunk number 2: Documento2_Oxysternonfestivum.Rnw:63-66
###################################################
plot(rbios,1,main="Cobertura arborea")
plot(Oshp,add=T,lwd=3)
 plot(CNEB,add=T,border="maroon")


###################################################
### code chunk number 3: Referencias
###################################################
if (!exists("scrb.cneb")) {
  load("~/NeoMapas/Rdata/MSM.rda")
  load("~/NeoMapas/Rdata/LIT.rda") ## hay una mejor version?
  (load("~/NeoMapas/Rdata/DatosLiteraturaScarabaeinae.rda"))
  scrb.cneb$genero <- sapply(scrb.cneb$especie,function(x) strsplit(x," ")[[1]][1])
  scrb.cneb$epiteto <- sapply(scrb.cneb$especie,function(x) strsplit(x," ")[[1]][2])
  scrb.cneb$lat[!is.na(scrb.cneb$lat) & scrb.cneb$lat < -5] <- scrb.cneb$lat[!is.na(scrb.cneb$lat) & scrb.cneb$lat < -5]*-1
  scrb.cneb$val <- paste(scrb.cneb$genero,scrb.cneb$epiteto)
}
cat(sprintf(" %1$s: \\citep{%1$s}; ",unique(scrb.cneb$cdg_ref)))


###################################################
### code chunk number 4: Documento2_Oxysternonfestivum.Rnw:97-104
###################################################
##SCB00033-00007 *O festivum
##SCB00033-00005 *O ebeninum
par(mar=c(0,0,0,0))
plot(vzla)
 points(lat~lon,subset(scrb.BDV,cdg_taxon %in% "SCB00033-00007"))
 points(lat~lon,subset(scrb.MIZA,grepl("Oxysternon f",scrb.MIZA$nombre)),col=2) 



###################################################
### code chunk number 5: Documento2_Oxysternonfestivum.Rnw:109-134
###################################################
 tmp1 <- Olit[,c("LON","LAT")]
  colnames(tmp1) <-c("lon","lat")
  
  tmp2 <- rbind(subset(scrb.BDV,cdg_taxon %in% "SCB00033-00007")[,c("lon","lat")],
                subset(scrb.MIZA,grepl("Oxysternon f",
                                       scrb.MIZA$nombre))[,c("lon","lat")])
  
  ##error en georef, punto en portuguesa realmente se refiere a el Pao en Edo Bolivar
  subset(scrb.BDV,lat>9 & grepl("festivum",nombre))
  subset(scrb.MIZA,lat>9 & grepl("festivum",nombre))
  
  ##
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat),"lon"] <- -62.6554584
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat),"lat"] <- 8.0345837
  
##rbind(tmp1,tmp2)
  Oll <- unique(rbind(tmp1,tmp2))
  Oll <- subset(Oll,!is.na(lat))
##plot(rbios,1)
##plot(vzla,add=T)

par(mar=c(0,0,0,0))
plot(vzla)
points(Oll)



###################################################
### code chunk number 6: Olt
###################################################
 tmp1 <- Olit[,c("LON","LAT")]
  colnames(tmp1) <-c("lon","lat")
  tmp1$pa <- TRUE
  
  tmp4 <- read.csv("~/NeoMapas/data/Eurysternus/EurysternusDB2012.csv",as.is=T,dec=",")[,c("Longitude","Latitude")]
  colnames(tmp4) <-c("lon","lat")
  tmp4$pa <- FALSE
  
  tmp2 <- scrb.BDV[,c("lon","lat")]
  tmp2$pa <- scrb.BDV$cdg_taxon %in% "SCB00033-00007"
  ##
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat) & tmp2$pa,"lon"] <- -62.6554584
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat) & tmp2$pa,"lat"] <- 8.0345837
  
  tmp3 <- scrb.MIZA[,c("lon","lat")]
  tmp3$pa <- grepl("Oxysternon f",scrb.MIZA$nombre)
  ##
  tmp3[tmp3$lat>9 & !is.na(tmp3$lat) & tmp3$pa,"lon"] <- -62.6554584
  tmp3[tmp3$lat>9 & !is.na(tmp3$lat) & tmp3$pa,"lat"] <- 8.0345837
  
  Olt <- rbind(tmp1,tmp2,tmp3,tmp4)
  ##dim(aggregate(Oll$pa,list(Oll$lon,Oll$lat),max))
  Olt <- unique(Olt)

  


###################################################
### code chunk number 7: Documento2_Oxysternonfestivum.Rnw:168-174
###################################################
par(mar=c(0,0,0,0))
Olt <- subset(Olt,!is.na(lon))
plot(rbios,1)
plot(vzla,add=T)
points(Olt[!Olt$pa,],pch=3,cex=.4)
points(Olt[Olt$pa,],pch=1,cex=.6,col=2)


###################################################
### code chunk number 8: coneccion
###################################################

if (!exists("trmp.NM")) {
 archivo.sql <- "~/NeoMapas/lib/sql/acceso.cnf"
 ##archivo.sql <- "~/NeoMapas/lib/sql/yapacana_restringido.cnf"
 
 mi.sql <- dbDriver("MySQL")
 
 
 c.nm2 <- dbConnect(mi.sql, group="Escarabajos",
                    default.file=archivo.sql)
 tmp1 <- dbGetQuery(c.nm2,"SET NAMES 'utf8'")
 tmp1 <- dbGetQuery(c.nm2,"SET CHARACTER SET utf8")
  
 prg <- paste("SELECT l.NM, YEAR(v.fecha1) as yr,",
              "DATE_FORMAT(t.fecha2,'%d/%m/%Y') as fch,",
              "v.cdg_trampa as vst,",
              "SUM((TIME_TO_SEC(TIMEDIFF(t.fecha2,",
              "t.fecha1)))/3600) as dur,",
              "Prd, lat, lon, alt ",
             "FROM TRAMPAS v",
              "LEFT JOIN LOCALIDADES l",
              "   on v.cdg_localidad=l.cdg_localidad",
              "  LEFT JOIN VISITAS t",
              "   on v.cdg_trampa=t.cdg_trampa",
              "  WHERE l.NM is not null",
             "  GROUP BY vst",
              "  ORDER BY vst, fch")
 trmp.NM <- dbGetQuery(c.nm2,prg)
 trmp.NM$sfrz <- trmp.NM$dur
 trmp.NM$sfrz.a <- unname(unlist(tapply(trmp.NM$sfrz,trmp.NM$vst,cumsum)))
 trmp.NM <- trmp.NM[trmp.NM$dur>0,]
 
  
 prg <- paste("SELECT l.NM, YEAR(v.fecha1) as yr,",
              "v.cdg_trampa as vst,",
              "t.cdg_tvn as tvn,",
              "cdg_ejemplar as jmp,",
              "e.cdg_taxonomia as spp,",
              "genero, especie, familia",
               "FROM EJEMPLARES j",
              "LEFT JOIN ESPECIE e",
              "   on e.cdg_taxon=j.cdg_especie",
              "  LEFT JOIN VISITAS t",
              "   on j.cdg_tvn=t.cdg_tvn",
             "LEFT JOIN TRAMPAS v",
              "   on v.cdg_trampa=t.cdg_trampa",
              "LEFT JOIN LOCALIDADES l",
              "   on v.cdg_localidad=l.cdg_localidad",
              "  WHERE familia like '%Scarabaeinae%'",
              "   && ((l.NM IN ('24','26','09','92')",
             "     && year(t.fecha1)<2006)",
              "     || (l.NM IN ('97') && year(t.fecha1)=2008)",
              "     || (l.NM IN ('93','05','08') && year(t.fecha1)=2006))")
 jmp.NM <- dbGetQuery(c.nm2,prg)
 
 prg <- paste("SELECT cdg_especie,l.NM,l.Prd,i.cdg_trampa,sum(n_ejemplares) as njmp,",
              "YEAR(t.fecha2) as yr",
              "FROM TABLA_ID i",
              "LEFT JOIN VISITAS t",
              "    ON t.cdg_tvn=i.cdg_tvn",
               "LEFT JOIN TRAMPAS v",
              "   on v.cdg_trampa=t.cdg_trampa",
              "LEFT JOIN LOCALIDADES l",
              "   on v.cdg_localidad=l.cdg_localidad",
               "GROUP BY cdg_especie,cdg_trampa")
 
 scrb.NM <- dbGetQuery(c.nm2,prg)
 scrb.NM$NMyr <- paste("NM",scrb.NM$NM," :: ",scrb.NM$yr,sep="")
 mtrz <- tapply(scrb.NM$njmp,list(scrb.NM$NMyr,scrb.NM$cdg_especie),sum)
  mtrz[is.na(mtrz)] <- 0
 dbDisconnect(c.nm2)
 
}
trmp.ll <- subset(trmp.NM,!is.na(lon) & lon<0 & lon >-75)

coordinates(trmp.ll) <- c("lon","lat")
proj4string(trmp.ll) <- CRS("+datum=WGS84 +proj=longlat")



###################################################
### code chunk number 9: DatosSolis
###################################################
load(file="~/NeoMapas/Rdata/scrbsolis.rda")

dtt <- trmp.ll$vst %in% scrb.solis$vst[scrb.solis$Oxysternon.festivum>0 | scrb.solis$Oxysternum.f..viridanum>0]
dt2 <- trmp.ll$vst %in% scrb.solis$vst[scrb.solis$Oxysternon.festivum>0]

tOxy <- aggregate(scrb.solis[,grep("Oxy",colnames(scrb.solis))],list(NM=sprintf("%02d",as.numeric(scrb.solis$NM))),sum)
rownames(tOxy) <- paste("NM",tOxy[,1],sep="")
tOxy <- tOxy[,-1]
colSums(tOxy)
##colSums(tOxy[,-1]>0)

fOxy <- aggregate(scrb.solis[,grep("Oxy",colnames(scrb.solis))],list(trampa=scrb.solis$vst),sum)
fOxy <- fOxy[fOxy$trampa %in% trmp.ll$vst,]

## Oxysternon para Arlene
##tOxy <- aggregate(scrb.solis[,grep("Oxy",colnames(scrb.solis))],list(trampa=scrb.solis$Grp),sum)
##colSums(tOxy[,-1])
##colSums(tOxy[,-1]>0)


###################################################
### code chunk number 10: Documento2_Oxysternonfestivum.Rnw:293-303
###################################################
##tOxy
##(trmp.ll@data[dtt,c("NM","vst")])

rowSums(tOxy)

##table(trmp.ll$NM[dtt])
##length(table(trmp.ll$vst[dtt]))
##dim(subset(fOxy,Oxysternon.festivum >0))
##subset(fOxy,Oxysternon.festivum >0)



###################################################
### code chunk number 11: Documento2_Oxysternonfestivum.Rnw:310-314
###################################################
par(mar=c(0,0,0,0))
plot(vzla,col="grey68",border="grey79")
points(trmp.ll,pch=3,cex=.3)
points(trmp.ll[dtt,],pch=1,cex=.4,col=2)


###################################################
### code chunk number 12: mediciones
###################################################
##setwd("/Users/jferrer/NeoMapas/etc/MedidasEscarabajos")
tmp <- read.csv("~/NeoMapas/etc/MedidasEscarabajos/20101014_MedicionesEscarabajosOFF.csv",as.is=T,dec=",")
##tmp$LC <- as.numeric(sub(",",".",tmp$LC))
##tmp$LE <- as.numeric(sub(",",".",tmp$LE))
##tmp$LF <- as.numeric(sub(",",".",tmp$LF))
##tmp$AP <- as.numeric(sub(",",".",tmp$AP))
tmp$fch <-  sapply(strsplit(tmp$Fecha,"/"),function(x){paste(x[3],x[2],x[1],sep="-")})
tmp$Sex <- factor(tmp$Sex)

##match(sprintf("NM%02d-T-%03d",tmp$NM,tmp$Trp),trmp.ll$vst)
tmpvrs <- extract(rbios,trmp.ll[match(sprintf("NM%02d-T-%03d",tmp$NM,tmp$Trp),trmp.ll$vst),])

##table(tmp$NM)



###################################################
### code chunk number 13: Documento2_Oxysternonfestivum.Rnw:339-341
###################################################
table(tmp$Trp)



###################################################
### code chunk number 14: SlcVars
###################################################
if (!exists("v0")) {
  v <- as.formula(paste("~", c("alt","arboreo",sprintf("bio%02d",1:19)), collapse="+"))
  v0 <- varclus(v, data=data.frame(values(sbios)))
}
plot(v0)
abline(h=.6,lty=2,col=2)
slc <- c("arboreo","bio01","bio02","bio03","bio05",
         "bio12", "bio18", "bio19")
slcbio <- rbios
values(slcbio) <- values(slcbio)[,slc]


###################################################
### code chunk number 15: mdl00
###################################################
if (!exists("me00")) {
   
  me00 <- maxent(slcbio,Oll)
  pr00 <- predict(me00,slcbio)
}


###################################################
### code chunk number 16: Documento2_Oxysternonfestivum.Rnw:382-385
###################################################
plot(pr00,col=brewer.pal(9,"BuGn"))

points(Oll)


###################################################
### code chunk number 17: mdl1
###################################################

if (!exists("me01")) {
  ## 
   prds <- data.frame(extract(rbios,Olt[,1:2]))
  
  me01 <- maxent(x=prds[,slc],p=Olt$pa+0)
  pr01 <- predict(me01,slcbio)
}



###################################################
### code chunk number 18: mapamejorado
###################################################
plot(pr01,col=brewer.pal(9,"BuGn"))
##plot(pr01>.5)
points(Olt,cex=0.5,pch=19,col=Oll$pa)
##plot(rbios,1)
plot(Oshp,add=T)
points(LAT~LON,Olit)

plot(vzla,add=T)
subset(Olt,lat>9 & pa)


###################################################
### code chunk number 19: test
###################################################
tst00 <- extract(pr00,trmp.ll)
tst01 <- extract(pr01,trmp.ll)
layout(matrix(1:2,ncol=2))
boxplot(tst00~dtt,main="modelo 1",ylab="Predicción",xlab="Observado")
boxplot(tst01~dtt,main="modelo 2")
##boxplot(extract(pr01,trmp.ll)~dt2)


###################################################
### code chunk number 20: AUC
###################################################
evaluate(tst00[dtt],tst00[!dtt])
evaluate(tst01[dtt],tst01[!dtt])


###################################################
### code chunk number 21: cortest
###################################################
cor.test(tst00,dtt+0)
cor.test(tst01,dtt+0)



###################################################
### code chunk number 22: otromapa
###################################################
plot(pr01,col=brewer.pal(9,"BuGn"))
points(subset(trmp.ll,!dtt),col=2,pch=3,cex=.5)
points(subset(trmp.ll,dtt),col=1,pch=19,cex=.5)


###################################################
### code chunk number 23: importancia
###################################################
plot(me01)
lines(vegan::bstick(me01, tot.var = 100,n=8),8:1,col=2)



###################################################
### code chunk number 24: response
###################################################
response(me01)



###################################################
### code chunk number 25: correlacion
###################################################

prd1 <- extract(pr01,trmp.ll)
vrs1 <- extract(rbios,trmp.ll)

## considerar solo las trampas que han sido revisadas>
##[match(fOxy$trampa,trmp.ll$vst),]



###################################################
### code chunk number 26: cordtt
###################################################
cor(vrs1[,slc],dtt+0,use="complete")
cor.test(vrs1[,slc[1]],dtt+0)


###################################################
### code chunk number 27: abundancia
###################################################
abnd <- fOxy$Oxysternon.festivum
abnd <- abnd[match(trmp.ll$vst,fOxy$trampa)]

cor(vrs1[,slc],log1p(abnd),use="complete")
cor.test(vrs1[,slc[1]],log1p(abnd))



###################################################
### code chunk number 28: Documento2_Oxysternonfestivum.Rnw:499-500
###################################################
cor.test(log1p(abnd),prd1,use="complete")


###################################################
### code chunk number 29: figcorrelacion (eval = FALSE)
###################################################
## ss <- extract(rlim,trmp.ll[match(fOxy$trampa,trmp.ll$vst),])
## ss <- prd1>.5
## plot(prd1[ss],log1p(abnd[ss]))
## ## limiter la prueba a la region donde se espera presencia (ver shapefile)
## ##cor.test(log1p(abnd)[!is.na(ss)],prd1[!is.na(ss)],use="complete")
## ##cor.test(log1p(abnd)[!is.na(ss)],boot::logit(prd1)[!is.na(ss)],use="complete")
## ##plot(boot::logit(prd1),log1p(abnd),subset=ss)
## ##cor.test(log1p(abnd)[ss],prd1[ss],use="complete")
## ##cor.test(log1p(abnd)[ss],vrs1[ss,"arboreo"],use="complete")
## ##cor.test(log1p(abnd)[ss],vrs1[ss,"arboreo"],use="complete")
## 
## ##plot(prd1,log1p(abnd),subset=ss)
## ##plot(prd1[ss],log1p(abnd[ss]))
## 
## ##write.csv(file="NeoMapas_Oxysternon_festivum.csv",data.frame(spp="Oxysternon festivum",trmp.ll@coords[match(unique(scrb.solis$vst[scrb.solis$Oxysternon.festivum>0 | scrb.solis$Oxysternum.f..viridanum>0]),trmp.xy@data$vst),]),row.names = FALSE)
## 
## ##write.csv(file="NeoMapas_Oxysternon_ebeninum.csv",data.frame(spp="Oxysternon ebeninum",trmp.ll@coords[match(unique(scrb.solis$vst[scrb.solis$Oxysternon.ebeninum>0]),trmp.xy@data$vst),]),row.names = FALSE)
## 
## ##write.csv(file="RevisionLiteratura_Oxysternon.csv",subset(scrb.BDV,grepl("Oxys",nombre)),row.names = FALSE)
## ##write.csv(file="RevisionMuseoMIZA_Oxysternon.csv",subset(scrb.MIZA,grepl("Oxys",nombre)),row.names = FALSE)
## 
## sfrz <- merge(aggregate(trmp.ll$sfrz,list(vst=trmp.ll$vst),sum,na.rm=T),
##          aggregate(trmp.ll@coords,list(vst=trmp.ll$vst),median,na.rm=T),by="vst")
## colnames(sfrz) <- c("codigo trampa","horas.trampa","lon","lat")
## 
## ##write.csv(file="NeoMapas_esfuerzoMuestreoScrb.csv",
## ##          sfrz,row.names = FALSE)
## 
## ##system("zip ScrbNeoMapasParaArlene.zip NeoMapas_esfuerzoMuestreoScrb.csv NeoMapas_Oxysternon_festivum.csv RevisionLiteratura_Oxysternon.csv NeoMapas_Oxysternon_ebeninum.csv RevisionMuseoMIZA_Oxysternon.csv")
## 


###################################################
### code chunk number 30: pares
###################################################
pairs(tmp[,c("LC","LE","LF","AP")])



###################################################
### code chunk number 31: Documento2_Oxysternonfestivum.Rnw:549-550
###################################################
print(xyplot(LE~AP|Sex,tmp))


###################################################
### code chunk number 32: Documento2_Oxysternonfestivum.Rnw:555-556
###################################################
print(xyplot(LC~AP,tmp))


###################################################
### code chunk number 33: Documento2_Oxysternonfestivum.Rnw:561-564
###################################################
##colSums(table(tmp$LC,tmp$Trp))
boxplot(tmp$LC~tmp$Trp,varwidth=T)



###################################################
### code chunk number 34: Documento2_Oxysternonfestivum.Rnw:568-570
###################################################
plot(tmp$LC~tmpvrs[,"arboreo"])



