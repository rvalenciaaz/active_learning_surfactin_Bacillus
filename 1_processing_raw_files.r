library(tidyverse)
library(rawrr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pracma)

folder_path<-"13_Oct_iter1_v2/"
f2<-str_sub(folder_path,end=-2)

names<-c()
for (i in  list.files(path = folder_path)){
    if (endsWith(i,".raw")){
        names<-append(names,i)
    }
}

extract<-function(name){
    rawfile<-paste(folder_path,name, sep="/")
    index<-rawrr::readIndex(rawfile = rawfile)
    
    scans<-unique(index$scanType)

    megatable <- data.frame(matrix(ncol = 3, nrow = 0))
    x <- c("rt", "tic","scan_type")
    colnames(megatable) <- x
    megatable$scan_type<-as.character(megatable$scan_type)

    for (i in scans){

      subi<-rawrr::readIndex(rawfile = rawfile) |>
        subset(scanType == i)
      listscan<-subi$scan
      S<- rawrr::readSpectrum(rawfile = rawfile, scan = listscan)
      for (j in 1:length(listscan)) {
        datos<-S[[j]]
        ticn<-datos$TIC
        sec<-datos$rtinseconds
        megatable<-megatable %>% add_row(rt = sec, tic = ticn,scan_type = i)
      }

    }
    megatable$file<-name
    return(megatable)
}

names

probar<-list()

for (j in 1:length(names)){
  probar[[j]] <- extract(names[j])
}

hypertable<-do.call("rbind", probar)

write_csv(hypertable, paste(f2,"_hypertable.csv",sep=""))



