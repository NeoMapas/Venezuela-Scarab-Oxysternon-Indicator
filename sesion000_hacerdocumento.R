##R --vanilla
if (system("hostname",intern=T) %in% c("tuitari","humboldt")) {
    options(java.parameters = "-Xmx11g" )
    setwd("~/tmp/NeoMapas")
    BIBDIR <- "/home/jferrer/Dropbox/NeoMapas/lib"
}
if (system("hostname",intern=T) %in% "catatumbo") {
    options(java.parameters = "-Xmx7g" )
    setwd("~/tmp/")
    BIBDIR <- "/home/mapoteca/Dropbox/NeoMapas/lib"
}

hoy <- format(Sys.time(), "%Y%m%d")
##install.packages("ROpenOffice", repos = "http://www.omegahat.org/R")
##install.packages(c("dismo","raster","ENMeval","AICcmodavg","pROC","compositions","tweedie"))
require(ROpenOffice)
require(dismo)
require(raster)
library(RMySQL)
require(gdata)
require(Hmisc)
require(lattice)
require(RColorBrewer)
require(ENMeval)
require(AICcmodavg)
require(MASS)
require(xtable)
require(pROC)
require(compositions)
require(vegan)
require(pscl)
require(plotrix)
require(tweedie)
library(statmod) # Needed to use  tweedie.profile

options(width=65)
options(width=65,stringsAsFactors=FALSE)


mi.path <- "~/NeoMapas/"
mi.path <- "~/Dropbox/NeoMapas/doc/"
mi.dir <- "902_AnalisisMorfometriaEscarabajos"
mi.arch <- "Documento1_MedicionesEscarabajos"
titulo <- "MedicionesEscarabajos"
mi.arch <- "Documento2_Oxysternonfestivum"
##read.csv("~/NeoMapas/etc/MedidasEscarabajos/20101021_MedicionesEscarabajosOFF.csv")
titulo <- "DistribucionAbundanciaBiomasas_del_TorococoVinotinto"


##mi.arch <- "Documento3_ScrbNM97"
##titulo <- "HojarasquitoAltosPipe"
source("~/Dropbox/Mapoteca/inc/inc00_funciones.R")
source(file=paste(mi.path,mi.dir,"/pre000_DataOfestivum.R",sep=""))


Onthophagus <- aggregate(scrb.solis[,grep("thophagus",colnames(scrb.solis))],by=list(vst=scrb.solis$vst),sum)

Trampas.e <- aggregate(data.frame(esfuerzo.muestreo=trmp.NM[,"sfrz"]),by=list(vst=trmp.NM$vst),sum)
Trampas.ll <- aggregate(trmp.NM[,c("lon","lat")],by=list(vst=trmp.NM$vst),mean)
Trampas <- merge(Trampas.e,Trampas.ll,by="vst")
Onthophagus <- merge(Trampas,Onthophagus,by="vst",all.x=T)

rownames(Onthophagus) <- Onthophagus$NM
Onthophagus <- Onthophagus[rowSums(Onthophagus[,-1])>0,-1]
Onthophagus <- subset(Onthophagus,lat>0 & lon< -40)

Onthophagus[is.na(Onthophagus)] <- 0

 elevacion <- raster("~/mapas/Venezuela/SRTM/SRTM.Venezuela.tif")
 bosque <- raster("~/mapas/Venezuela/VCF/MOD44B.2006.Venezuela.TREE.tif")
Onthophagus$bosque <- extract(bosque,Onthophagus[,c("lon","lat")])
Onthophagus$elevacion <- extract(elevacion,Onthophagus[,c("lon","lat")])
Onthophagus$elevacion[Onthophagus$elevacion<0] <- 0
Onthophagus$bosque[Onthophagus$bosque>100] <- NA

save(file="~/Dropbox/NeoMapas/Rdata/EjercicioOnthophagus.rda",Onthophagus)

## ver ejemplo Area de Referencia: AreaReferencia_Ofestivum.R

## primero hay que correr Documento2_Oxysternonfestivum
mi.arch <- "Documento2_Oxysternonfestivum"
mi.arch <- "Documento20_OnlineResourceJICO"

##system(sprintf("rm %s*",mi.arch))
system(sprintf("sed s:BIBDIR:%s:g %s%s/%s.Rnw > %s.Rnw",BIBDIR,mi.path,mi.dir,mi.arch,mi.arch))
##Sweave(file=paste(mi.path,mi.dir,"/",mi.arch,".Rnw",sep=""),eps=F)
Sweave(file=sprintf("%s.Rnw",mi.arch),eps=F)
##Stangle(file=paste(mi.path,mi.dir,"/",mi.arch,".Rnw",sep=""))
Stangle(file=sprintf("%s.Rnw",mi.arch))
tools::texi2dvi(paste(mi.arch,".tex",sep=""), pdf=TRUE)

ssystem(sprintf("evince %s.pdf &",mi.arch))

system(paste("mv ",mi.arch,".pdf ",mi.path,"/",mi.dir,"/",hoy,"_",titulo,".pdf",sep=""))
system(paste("mv ",mi.arch,".R ",mi.path,"/",mi.dir,"/",hoy,"_",titulo,".R",sep=""))



#####################
### lista de escarabajos por transecta
#############
  mi.arch <- "ListaEspeciesTodas"
  system(sprintf("rm %s",mi.arch))

for (nn in as.numeric(unique(scrb.solis$NM))) {
  ss <- subset(scrb.solis,NM %in% nn)

  tt <- with(ss,aggregate(ss[,6:158],
                          list(fecha=yr,grupo=Grp,trampa=vst),
                          sum))
  tt <- tt[order(tt$fecha,tt$grupo,tt$trampa),c(T,T,T,colSums(tt[,-(1:3)])>0)]
  
##  head(tt)

  dts <- data.frame(fecha=rep(tt$fecha,ncol(tt)-3),Grupo=rep(tt$grupo,ncol(tt)-3),trampa=rep(tt$trampa,ncol(tt)-3),stack(tt[,-(1:3)]))
  dts <- subset(dts,values>0)

##  mi.arch <- sprintf("ListaEspeciesNM%02i",nn)
  cat(sprintf("
::::::::.::::::::::.:::::::::::::::::.::::::::::::::::::.:::::::::
Lista de especies de escarabajos para la transecta NM%02i
Según identificación realizada por Ángel Solís, INBio, Costa Rica
::::::::.::::::::::.:::::::::::::::::.::::::::::::::::::.:::::::::\n",nn),
      file=mi.arch,append=T)
  for (j in levels(dts$ind)) {
    ss <- subset(dts,ind %in% j)
    cat(sprintf("\n%s",j),file=mi.arch,append=T)
    for (y in unique(ss$fecha)) {
      cat(sprintf("\nAño %s:",y),file=mi.arch,append=T)
      for (g in unique(ss$Grupo)) {
        kk <- subset(ss,fecha %in% y & Grupo %in% g)
        if (nrow(kk)>0)
          cat(sprintf(" (Grupo %s) trampas %s.",
                      g,paste(kk$trampa,collapse=", ")),file=mi.arch,append=T)
      }
    }
    cat(sprintf("\n"),file=mi.arch,append=T)
    
  }
}

##lr.stat <- lrtest(mod0, mod1)
##(1-exp(-as.numeric(lr.stat$stats[1])/n))/(1-exp(2*as.numeric(logLik(mod0)/n)))


unlist(lapply(PS.lst,function(x) summary(x)$adj.r.squared))
unlist(lapply(AP.lst,function(x) summary(x)$adj.r.squared))
unlist(lapply(LE.lst,function(x) summary(x)$adj.r.squared))
unlist(lapply(LF.lst,function(x) summary(x)$adj.r.squared))


##n <- nrow(lms[[1]]$model)
##null.dev <- -2*logLik(lms[["Psi(HS)"]])
##mi.dev <- unlist(lapply(lms[cvg],function(x) -2*logLik(x)))
##R2 <- (1 - exp((mi.dev - null.dev)/n))/(1 - exp(-null.dev/n))


##n <- PS.lst[["PS(.)"]]$df.residual
##L0 <- logLik(PS.lst[["PS(.)"]])
##LR <- unlist(lapply(PS.lst,function(x) 2*(logLik(x)-L0)))
##(1-exp(-LR/n))/(1-exp(-(-2*L0)/n))


##Oxysternon en Venezuela
##festivum
##conspiciliatum
##ebeninum
##silenus (= smargdinum)
##spiniferum

table(subset(scrb.MIZA,grepl("Oxys",nombre))$nombre)
table(subset(scrb.BDV,grepl("Oxys",nombre))$nombre)
with(subset(scrb.BDV,grepl("Oxys",nombre)),table(cdg_ref,nombre))

subset(fOxy,Oxysternum.f..viridanum>0)
subset(fOxy,Oxysternon.ebeninum>0)
