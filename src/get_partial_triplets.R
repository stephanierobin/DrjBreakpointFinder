# #############################################################################
# This file contains the definition of the function getPartialTriplets
#
# Author: Claire Lemaitre
# #############################################################################


## 28/08/2012
source("hits_functions.R")

## exemple :
##  getPartialTriplets("../res-sanger-full/blast/megablast-default.m8","../data/formatted-data/HdIV_SANGER_clean_min100_rlength.tab","../res-sanger-full/drjPairs_confirmed.tab","test/seqCoord.tab","test/drjPair2id.tab","../res-sanger-full/breakpoint/seqCoordPairs.tab")


## Fonction pour obtenir les triplets : Read-DRJ1-DRJ2 des reads qui chevauchent seulement partiellement des DRJs déjà connues (il suffit qu'une seule des 2 DRJs soit couverte par un hit avec ce read)
## Paramètres :
##  - blastBrutFile : fichier brut de sortie de blast ou megablast (format m8)
##  - readLengthFile : fichier contenant la taille de chacun des reads
##  - bestReadsFile : fichier contenant les reads déjà utilisés, par exemple selectpairs.tab
##  - DRJpairsFile : fichier avec les coord des paires de DRJs connues
##  - minPcId : %id mini du hit pour être retenu
##  - minDRJcov : couverture minimale de la DRJ pour être retenu
## Renvoie un tableau avec les triplets : id bac read inf1 sup1 inf2 sup2 orient (bornes des DRJs sur le bac)
getPartialTriplets=function(blastBrutFile,readLengthFile,DRJpairsFile,coordFile,drjPair2IdsFile,bestReadsFile=NULL,id2drjFile=NULL,minPcId=90,minDRJcov=50,minReadCov=85,extend=0){

  blast=formatBlastTable(blastBrutFile,readLengthFile)

  ## remove_best_reads
  ## enlève les reads pris dans la methode 1 = reads qui chevauchent entierement les DRJs
  if(!is.null(bestReadsFile)){
    bestReads=as.character(read.table(bestReadsFile,h=T)$read)
    blast=blast[!is.element(blast$read,bestReads),]
  }

  ## enleve les hits à faible %id :
  blast=blast[blast$pcid>=minPcId,]

  ## enleve les hits qui couvrent trop peu le read
  blast=blast[blast$length>=minReadCov*blast$l/100,]

  if(nrow(blast)>0){
    ## selectionne les paires de DRJs pour chaque hit de chaque read
    DRJpairs=read.table(DRJpairsFile,h=T,fill=T)
    triplets=apply(blast[,c("read","bac","inf2","sup2","orient")],1,function(x) getDRJpairs(x,DRJpairs,minDRJcov))
    triplets=unique(do.call("rbind",triplets))

    triplets$id=1:nrow(triplets)
    triplets=triplets[,c("id","read","bac","inf1","sup1","inf2","sup2","orient")]

    ## id2drj
    if(!is.null(id2drjFile)){
      write.table(triplets[,c("id","bac","inf1","sup1","inf2","sup2")],id2drjFile,quote=F,row.names=F)
    }

    ## obtenir les idList pour chaque paire de DRJ :
    res=by(triplets,paste(triplets$bac,triplets$inf1,triplets$sup1,triplets$inf2,triplets$sup2,sep="-"), function(x) {
      nbHits=nrow(x)
      nbReads=length(unique(x$read))
      idList=paste(x$id,collapse=",")
      return(as.character(c(as.character(x$bac[1]),x$inf1[1],x$sup1[1],x$inf2[1],x$sup2[1],nbHits,nbReads,idList)))
    })

    mat=matrix(unlist(res),ncol=8,byrow=T)
    final=as.data.frame(mat)
    names(final)=c("bac","inf1","sup1","inf2","sup2","nbHits","nbReads","idList")

    write.table(final,drjPair2IdsFile,quote=F,row.names=F)

    ## obtenir les coordonnées des séquences sur le bac, soit les coord des drjs, soit on peut les étendre à gauche et à droite de la valeur "extend"
    coord=triplets
    coord$inf1=pmax(coord$inf1-extend,1)
    coord$sup1=coord$sup1+extend
    coord$inf2=pmax(coord$inf2-extend,1)
    coord$sup2=coord$sup2+extend

    write.table(coord,coordFile,quote=F,row.names=F)
  }
  else{
    print("No significative hits found on known DRJs")
  }
}



## Pour un hit de blast récupère les paires de DRJs pour lesquelles le hit chevauche d'au moins minDRJcov % au moins une des 2 DRJs
getDRJpairs=function(blastLine,DRJpairs,minDRJcov){

  bac=blastLine[2]
  infR=as.numeric(blastLine[3])
  supR=as.numeric(blastLine[4])

  tab=DRJpairs[as.character(DRJpairs$bac)==bac,]
  if(nrow(tab)>0){
    tab$minT=(tab$sup1-tab$inf1)*minDRJcov/100

    ## over.sup - over.inf : si positif = taille de l'overlap (si négatif = pas d'overlap)
    over1.inf=pmax(tab$inf1,infR)
    over1.sup=pmin(tab$sup1,supR)
    over2.inf=pmax(tab$inf2,infR)
    over2.sup=pmin(tab$sup2,supR)

    select=tab[(over1.sup-over1.inf+1)>=tab$minT | (over2.sup-over2.inf+1)>=tab$minT,]

    select$read=rep(blastLine[1],nrow(select))
    select$orient=rep(blastLine[5],nrow(select))
    return(select)
  }
}
