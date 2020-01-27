#################
### SIG
################
###
## Variables ambientales para modelo dist.
###
if (!exists("rbios")) {
    TMWB <- shapefile("~/gisdata/vectorial/TMworldborders/TM_WORLD_BORDERS-0.3.shp")
    
    load("~/Dropbox/NeoMapas/Rdata/SIG.rda")
    rbios <- stack(c("~/Dropbox/NeoMapas/data/gisdata/Ofestivum/alt.tif",sprintf("~/Dropbox/NeoMapas/data/gisdata/Ofestivum/bio%s.tif",1:19)))
    values(rbios)[values(rbios)>55000] <- NA
    rarbo <- raster(sprintf("~/Dropbox/NeoMapas/data/gisdata/Ofestivum/arboreo.tif"))
  values(rarbo)[values(rarbo)==0] <- NA
  values(rarbo)[values(rarbo)==254] <- 0
  values(rarbo)[values(rarbo)==255] <- 5
  r1 <- aggregate(rarbo,5,fun=mean,na.rm=TRUE)
  arboreo <- raster::resample(r1, raster(rbios,1), method="bilinear")
  rbios <- stack(arboreo,rbios)
  Oshp <- shapefile("~/Dropbox/NeoMapas/data/Oxysternon/Ofestivum.shp")
  ##bbox -73.3755402 -46.30876 -3.470826 13
  rlim <- rasterize(Oshp,arboreo)
  sbios <- rbios*rlim
  names(sbios) <- c("arboreo","alt",sprintf("bio%02d",1:19))
  rbios <- stack(rbios,!is.na(rlim))
  names(rbios) <- c("arboreo","alt",sprintf("bio%02d",1:19),"M")
}





################################
## Datos de distribución
#################################



##<<GBIF>>=
require(rgbif)

## 17.000 registros?
##scrb.gb <- read.table("~/Dropbox/NeoMapas/data/Oxysternon/0000294-150827100048397.csv",sep="\t",header=T)

##scrb.gb <- read.ods("~/Dropbox/NeoMapas/data/Oxysternon/0000294-150827100048397.ods")
##save(scrb.gb,file="~/Dropbox/NeoMapas/Rdata/ScarabaeidaeGBIF.rda")
(load("~/Dropbox/NeoMapas/Rdata/ScarabaeidaeGBIF.rda"))
scrb.gb$lat <- as.numeric(unlist(scrb.gb$decimallatitude))
scrb.gb$lon <- as.numeric(unlist(scrb.gb$decimallongitude))
scrb.gb$lon[scrb.gb$lon < -100] <- scrb.gb$lon[scrb.gb$lon < -100]/1000
scrb.gb$lat[scrb.gb$lat < -100] <- scrb.gb$lat[scrb.gb$lat < -100]/1000
scrb.gb$lat[scrb.gb$lat > 100] <- scrb.gb$lat[scrb.gb$lat > 100]/1000

coordinates(scrb.gb) <- c("lon","lat")
proj4string(scrb.gb) <- CRS("+datum=WGS84 +proj=longlat")

### intento con rgbif tiene problemas para bajar todos los registros
if (!exists("scrb.gb")) {
    qry <- occ_search(scientificName="Oxysternon",limit=20000)
    oxyf.gb <- subset(qry$data,!is.na(decimalLongitude))
    coordinates(oxyf.gb) <- c("decimalLongitude","decimalLatitude")
    proj4string(oxyf.gb) <- CRS("+datum=WGS84 +proj=longlat")


    qry <- occ_search(scientificName="Scarabaeinae",limit=20000)
    e <- extent(rbios)
    qr2 <- occ_search(taxonKey="5840",geometry=sprintf('POLYGON((%1$s %3$s, %1$s %4$s, %2$s %4$s, %2$s %3$s, %1$s %3$s))',e[1], e[2],e[3],e[4]), limit=200000)
 
}


table(subset(scrb.gb@data,grepl("festivum",scientificname))$scientificname)

###
## Referencias bibliográficas
###
if (!exists("scrb.cneb")) {
  load("~/Dropbox/NeoMapas/Rdata/MSM.rda")
  load("~/Dropbox/NeoMapas/Rdata/LIT.rda") ## hay una mejor version?
  (load("~/Dropbox/NeoMapas/Rdata/DatosLiteraturaScarabaeinae.rda"))
  scrb.cneb$genero <- sapply(scrb.cneb$especie,function(x) strsplit(x," ")[[1]][1])
  scrb.cneb$epiteto <- sapply(scrb.cneb$especie,function(x) strsplit(x," ")[[1]][2])
  scrb.cneb$lat[!is.na(scrb.cneb$lat) & scrb.cneb$lat < -5] <- scrb.cneb$lat[!is.na(scrb.cneb$lat) & scrb.cneb$lat < -5]*-1
  scrb.cneb$val <- paste(scrb.cneb$genero,scrb.cneb$epiteto)
}


###
## Edmonds
###
  Olit <- read.ods(file="~/Dropbox/NeoMapas/data/Oxysternon/LocalidadesOfestivum.ods")

####
## Revisión MALUZ
###

##07°11'35''N - 05°00'08" ??
## ver http://humbertosilvacubillan.blogspot.com/2010/04/el-caura-reserva-forestal-en-la.html
##7.342567	-65.211855	301	Unavailable
##6.003105	-61.392155	3036	Unavailable
scrb.MALUZ <- data.frame(lat=c(7+(11+(35/60))/60,6.003105),
lon=c(-1*(65+(0+(8/60))/60),-61.392155),
especie="Oxysternon festivum")


###
## resumen datos Vzla

## otra lit
tmp1 <- Olit[,c("LON","LAT")]
colnames(tmp1) <-c("lon","lat")

##literatura y colecciones
tmp2 <- rbind(subset(scrb.BDV,cdg_taxon %in% "SCB00033-00007")[,c("lon","lat")],
              subset(scrb.MIZA,grepl("Oxysternon f",
                                     scrb.MIZA$nombre))[,c("lon","lat")],
              scrb.MALUZ[,c("lon","lat")])

tmp3 <- coordinates(scrb.gb)[grepl("festivum",scrb.gb@data$scientificName),]
colnames(tmp3) <- c("lon","lat")
##error en georef, punto en portuguesa realmente se refiere a el Pao en Edo Bolivar
  subset(scrb.BDV,lat>9 & grepl("festivum",nombre))
  subset(scrb.MIZA,lat>9 & grepl("festivum",nombre))
  
  ##
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat),"lon"] <- -62.6554584
  tmp2[tmp2$lat>9 & !is.na(tmp2$lat),"lat"] <- 8.0345837
  
##rbind(tmp1,tmp2)
  Oll <- unique(rbind(tmp1,tmp2,tmp3))
  Oll <- subset(Oll,!is.na(lat))
##plot(rbios,1)
##plot(vzla,add=T)

 tmp1 <- Olit[,c("LON","LAT")]
  colnames(tmp1) <-c("lon","lat")
  tmp1$pa <- TRUE
  
  tmp4 <- read.csv("~/Dropbox/NeoMapas/data/Eurysternus/EurysternusDB2012.csv",as.is=T,dec=",")[,c("Longitude","Latitude")]
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

tmp5 <- data.frame(coordinates(scrb.gb))
colnames(tmp5) <- c("lon","lat")
tmp5$pa <- grepl("festivum",scrb.gb@data$scientificname)

  Olt <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5)
  ##dim(aggregate(Oll$pa,list(Oll$lon,Oll$lat),max))
  Olt <- unique(Olt)

###
## datos NeoMapas (DB)
###
(load("~/Dropbox/NeoMapas/Rdata/NMscrb.rda"))

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
trmp.NM$year <- as.numeric(substr(trmp.NM[,"fch"],1,4)) 
trmp.NM$month <- as.numeric(substr(trmp.NM[,"fch"],6,7)) 
trmp.NM$doy <- cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31))[as.numeric(substr(trmp.NM[,"fch"],6,7))-1] + 
        as.numeric(substr(trmp.NM[,"fch"],9,10))

## error detectado por valores atípicos de LAI, EVI, etc, corregimos:
trmp.NM[trmp.NM$vst %in% "NM39-T-034","lat"] <- 8.52652

if (!exists("trmp.ll")) {
    trmp.ll <- subset(trmp.NM,!is.na(lon) & lon<0 & lon >-75)
    trmp.ll <- merge(merge(aggregate(trmp.ll[,c("njmp","sfrz")],
                                     list(NM=trmp.ll$NM,Prd=trmp.ll$Prd,
                                          vst=trmp.ll$vst),sum),
                           aggregate(trmp.ll[,c("year","month","doy")],
                                     list(NM=trmp.ll$NM,Prd=trmp.ll$Prd,
                                          vst=trmp.ll$vst),function(x) round(mean(x))),
                           by=c("NM","Prd","vst")),
                     aggregate(trmp.ll[,c("lon","lat")],
                               list(NM=trmp.ll$NM,Prd=trmp.ll$Prd,
                                    vst=trmp.ll$vst),mean),by=c("NM","Prd","vst"))
    
    trmp.ll$fPET <- sprintf("%sM%02d",trmp.ll$year,trmp.ll$month)
    trmp.ll$fLC1 <- sprintf("A%s",trmp.ll$year)
    trmp.ll$fEVI <- sprintf("%s%s",trmp.ll$year,as.character(cut(trmp.ll$doy,breaks= c(seq(1,365,by=16),365),label=sprintf("%03d",seq(1,365,by=16)))))
    trmp.ll$fLST <- sprintf("%s%s",trmp.ll$year,as.character(cut(trmp.ll$doy,breaks= c(seq(1,365,by=8),365),label=sprintf("%03d",seq(1,365,by=8)))))
    
    coordinates(trmp.ll) <- c("lon","lat")
    proj4string(trmp.ll) <- CRS("+datum=WGS84 +proj=longlat")
}

## variables MODIS

rslt <- trmp.ll@data
xys <- data.frame(trmp.ll@coords)
fiX <- "NM"
fiY <- c("year","month","doy")
source("~/Dropbox/Mapoteca/inc/inc03_extractModis.R",echo=T)



for (vv in c("LST_Day_1km","LST_Night_1km","PET_1km","ET_1km",
             "Lai_1km","Fpar_1km","v250m_16_days_NDVI","v250m_16_days_EVI")) 
    trmp.ll@data[,vv] <- na.fix[,vv]

##stop("ahora el QC")
boxplot(trmp.ll@data$Fpar_1km ~ rslt$FparLai_QC)

## chequeo
colSums(is.na(trmp.ll@data))
trmp.ll@data[trmp.ll@data$vst %in% c("NM66-TRP070","NM66-TRP071"),"LST_Night_1km"] <- rep(trmp.ll@data[trmp.ll@data$vst %in% "NM66-TRP069","LST_Night_1km"],2)

trmp.ll@data$DET <- trmp.ll@data$ET_1km/trmp.ll@data$PET_1km


##error en coordenadas de NM39-T-034 resultaba en valores muy bajos de EVI, Lai, etc, corregimos arriba y chequeamos aquí:

trmp.ll@data[ss & trmp.ll@data$v250m_16_days_EVI <0,]
subset(trmp.ll@data,vst %in% c("NM39-T-034","NM39-T-035","NM39-T-033"))
subset(coordinates(trmp.ll),trmp.ll@data$vst %in% c("NM39-T-034","NM39-T-035","NM39-T-033"))


##scrb.solis
load(file="~/Dropbox/NeoMapas/Rdata/scrbsolis.rda")
## corregir nombre de trampa de Interceptación de vuelo
scrb.solis[scrb.solis$Prd %in% "TI-2" & scrb.solis$NM %in% "41","vst"] <- "NM41-I-002"

## O. festivum identificados por Solis
##dtt <- trmp.ll$vst %in% scrb.solis$vst[scrb.solis$Oxysternon.festivum>0 | scrb.solis$Oxysternum.f..viridanum>0]
##dtt <- trmp.ll$vst %in% scrb.solis$vst[scrb.solis$Oxysternon.festivum>0]



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

#########
## <<mediciones,echo=false>>=
###################
##setwd("/Users/jferrer/NeoMapas/etc/MedidasEscarabajos")
##mds <- read.csv("~/NeoMapas/etc/MedidasEscarabajos/20150225_MedicionesEscarabajosOFF.csv",as.is=T,dec=",")
##mds <- read.csv("~/Dropbox/NeoMapas/data/MedidasEscarabajos/20150421_MedicionesEscarabajosOFF.csv",as.is=T,dec=",",sep="\t")
##mds <- read.csv("~/Dropbox/NeoMapas/data/MedidasEscarabajos/20150828_MedicionesEscarabajosOFF.csv",as.is=T,dec=",",sep=",")
mds <- read.csv("~/Dropbox/NeoMapas/data/MedidasEscarabajos/20150925_MedicionesEscarabajosOFF.csv",as.is=T,dec=",",sep=",")
##mds <- read.csv("~/NeoMapas/etc/MedidasEscarabajos/20101014_MedicionesEscarabajosOFF.csv",as.is=T,dec=",")
##mds$LC <- as.numeric(sub(",",".",mds$LC))
##mds$LE <- as.numeric(sub(",",".",mds$LE))
##mds$LF <- as.numeric(sub(",",".",mds$LF))
##mds$AP <- as.numeric(sub(",",".",mds$AP))
mds$fch <-  sapply(strsplit(mds$Fecha,"/"),
                   function(x){ifelse(length(x)==3,sprintf("%s-%02i-%02i",
                                       x[3],as.numeric(x[2]),
                                       as.numeric(x[1])),NA)})
mds$Sex <- factor(mds$Sex)

## colocamos vst directamente en el archivo (ver. 2015-09-23)
##mds$vst <- sprintf("NM%02d-T-%03d",as.numeric(mds$NM),as.numeric(mds$Trp))
##ss <- mds$NM %in% c(24) & grepl(2009,mds$Fecha)
##mds$vst[ss] <- sprintf("NM%02d-TR-%03d",
##                       as.numeric(mds$NM[ss]),
##                       as.numeric(mds$Trp[ss]))
##ss <- mds$NM %in% c(22,93) & grepl(2006,mds$Fecha)
##mds$vst[ss] <- sprintf("%02d-T-%03d",
##                       as.numeric(mds$NM[ss]),
##                       as.numeric(mds$Trp[ss]))
##ss <- mds$NM %in% c(34,65) 
##mds$vst[ss] <- sprintf("TRP%02dt%03d",
##                       as.numeric(mds$NM[ss]),
##                       as.numeric(mds$Trp[ss]))
##ss <- mds$Trp %in% "TIV-2"
##mds$vst[ss] <- "NM41-I-002"
##subset(mds,!vst %in% trmp.ll$vst)
##mds <- subset(mds,vst %in% trmp.ll$vst)


## adiciones de ejemplares identificados por CL al momento de las mediciones
## los de NM41 no fueron revisados por Solis en su última visita
## chequear "NM39-T-002" "NM39-T-035"

## la lista revisada tiene todos los ejemplares medidos y aquellos reportados por Solis y Ascencao y aún no medidos
dtt <- trmp.ll$vst %in% mds$vst

## adiciones de ejemplares en base de datos pero no en colección
## muestreo de 2005 en NM24 y trampa de interceptación de vuelo en NM39
## ahora reflejado en el archivo de mediciones (ver. 2015-09-23)
##dtt[trmp.ll$vst %in% subset(scr.NM,genero %in% "Oxysternon" & especie %in% "festivum")$vst] <- T

##match(sprintf("NM%02d-T-%03d",mds$NM,mds$Trp),trmp.ll$vst)
mdsvrs <- extract(rbios,trmp.ll[match(mds$vst,trmp.ll$vst),])

##mdsvr2 <- trmp.ll@data[match(mds$vst,trmp.ll$vst),c("LST_Day_1km","PET_1km","LST_Night_1km","ET_1km","Lai_1km","Fpar_1km","v250m_16_days_EVI","v250m_16_days_NDVI","HS")]
mdsvr2 <- trmp.ll@data[match(mds$vst,trmp.ll$vst),c("LST_Day_1km","PET_1km","LST_Night_1km","ET_1km","Lai_1km","Fpar_1km","v250m_16_days_EVI","v250m_16_days_NDVI","DET")]


##table(mds$NM)
mds$sgrps <- NA
mds$sgrps[mds$Sex %in% "♀"] <- "H"
mds$sgrps[mds$Sex %in% "♂"] <- "M1"
mds$sgrps[mds$Sex %in% "♂" & mds$LC < 0.4] <- "M2"


##xyplot(PS~LE|sgrps,mds)




##write.csv(file="NeoMapas_Oxysternon_festivum.csv",data.frame(spp="Oxysternon festivum",trmp.ll@coords[match(unique(scrb.solis$vst[scrb.solis$Oxysternon.festivum>0 | scrb.solis$Oxysternum.f..viridanum>0]),trmp.xy@data$vst),]),row.names = FALSE)

##write.csv(file="NeoMapas_Oxysternon_ebeninum.csv",data.frame(spp="Oxysternon ebeninum",trmp.ll@coords[match(unique(scrb.solis$vst[scrb.solis$Oxysternon.ebeninum>0]),trmp.xy@data$vst),]),row.names = FALSE)

##write.csv(file="RevisionLiteratura_Oxysternon.csv",subset(scrb.BDV,grepl("Oxys",nombre)),row.names = FALSE)
##write.csv(file="RevisionMuseoMIZA_Oxysternon.csv",subset(scrb.MIZA,grepl("Oxys",nombre)),row.names = FALSE)

sfrz <- merge(aggregate(trmp.ll$sfrz,list(vst=trmp.ll$vst),sum,na.rm=T),
         aggregate(trmp.ll@coords,list(vst=trmp.ll$vst),median,na.rm=T),by="vst")
colnames(sfrz) <- c("codigo trampa","horas.trampa","lon","lat")

##write.csv(file="NeoMapas_esfuerzoMuestreoScrb.csv",
##          sfrz,row.names = FALSE)

##system("zip ScrbNeoMapasParaArlene.zip NeoMapas_esfuerzoMuestreoScrb.csv NeoMapas_Oxysternon_festivum.csv RevisionLiteratura_Oxysternon.csv NeoMapas_Oxysternon_ebeninum.csv RevisionMuseoMIZA_Oxysternon.csv")




###########
## Para la presentación
##############

##png("~/Escritorio/PresentacionOxysternon_MapaDatos.png",width=700,height=500)
##plot(rlim,legend=F)
##plot(vzla,add=T)
##points(Olt,pch=3,cex=.3)
##points(subset(Olt,pa),pch=19,col=2,cex=.8)
##dev.off()

##png("~/Escritorio/PresentacionOxysternon_MapaNeoMapas.png",width=700,height=500)
##plot(rlim,legend=F)
##plot(vzla,add=T)
##points(trmp.ll,pch=3,cex=.3)
##points(trmp.ll[dtt,],pch=19,col=2,cex=.8)
##dev.off()

##png("~/Escritorio/PresentacionOxysternon_MaxEntRandomBackground.png",width=700,height=500)
##plot(prA13,legend=T,breaks=round(seq(0,1,length=10),2),col=brewer.pal(9,"Spectral"))
##plot(vzla,add=T,border="yellow")
##plot(Oshp,add=T)
##points(Olt,pch=3,cex=.3)
##points(subset(Olt,pa),pch=19,col=2,cex=.8)
##dev.off()


##png("~/Escritorio/PresentacionOxysternon_MaxEntSampleBackground.png",width=700,height=500)
##plot(prJ30,legend=T,breaks=round(seq(0,1,length=10),2),col=brewer.pal(9,"Spectral"))
##plot(vzla,add=T,border="yellow")
##plot(Oshp,add=T)
##points(Olt,pch=3,cex=.3)
##points(subset(Olt,pa),pch=19,col=2,cex=.8)
##dev.off()




##########
## cuales faltan por revisar

unq <- unique(trmp.NM[,c("NM","yr")])
unw <- unique(scrb.solis[,c("NM","yr")])
unw$NM <- sprintf("%02d",as.numeric(unw$NM))
unw$yr <- sprintf("20%s",unw$yr)

apply(unq,1,paste,collapse="::")
unq[!apply(unq,1,paste,collapse="::") %in% apply(unw,1,paste,collapse="::"),]

x2 <- unique(subset(scrb.solis,Oxysternon.festivum>0 & NM %in% 41)[,c("NM","vst","Prd")])
x1 <- unique(subset(mds,NM %in% 41)[,c("NM","Trp","vst")])
merge(x1,x2,by=c("NM","vst"),all=T)

