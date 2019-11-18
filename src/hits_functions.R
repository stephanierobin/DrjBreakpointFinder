## Claire Lemaitre
## 03/08/2011

## some functions for parsing blast hits


## Formate le fichier brut de résultat de blast/megablast (m8)
##  - enleve les colonnes inutiles
##  - rajoute la colonne "orient" et ramène tous les inf<sup
##  - rajoute la colonne "l" : taille du read (utilise le fichier readLengthFile
formatBlastTable=function(blastBrutFile,readLengthFile=NULL){

  tab=read.table(blastBrutFile)
  tab=tab[,-c(5,6,12)]## supprime les colonnes mistmatches, gap openings et bit score
  names(tab)=c("read","bac","pcid","length","inf1","sup1","inf2","sup2","evalue")
  tab$couple=paste(tab$read,tab$bac,sep="---")
  tab$orient=ifelse(tab$inf2<tab$sup2,1,-1)

  ## met toujours inf2<sup2
  inf2=pmin(tab$inf2,tab$sup2)
  sup2=pmax(tab$inf2,tab$sup2)
  tab$inf2=inf2
  tab$sup2=sup2

  ## rajoute la colonne l : taille du read
  if(!is.null(readLengthFile)){
    lengthTab=read.table(readLengthFile,h=FALSE,colClasses=c("character","integer"))
    names(lengthTab)=c("read","l")
    tab=merge(tab,lengthTab,by="read",sort=F)
  }

  tab$id=1:nrow(tab)
  tab=tab[,c(ncol(tab),1:(ncol(tab)-1))]

  return(tab)
}


## removing embedded hits :

## we consider a table with only forward hits (orient=1) inf<sup for both query and target
removeEmbeddedHits=function(blastTable,self=F,rename=F){
  ## embedded hit means that :
  ## hit 1   ----------            --------
  ## hit 2     ------              -----
  ## both segments on query and target are embedded in both segments of the first hit

  ## options :
  ##  -self =T : query=target
  ##  -rename =T : returns the filtered table with names inf1,sup1, etc. WARNING : in this table inf<sup (even for reverse hits)

  if(!is.element("couple",names(blastTable))){
    if(self){
      couples=blastTable$V1
    }
    else{
      couples=paste(blastTable$V1,blastTable$V2,sep="---")
    }
  }
  else{
    couples=blastTable$couple
  }
  
  if(!is.element("id",names(blastTable))){
    blastTable$id=1:nrow(blastTable)
  }

  tmpTable=blastTable
  ## rename des colonnes de position
  if(!is.element("inf1",names(blastTable))){
    names(tmpTable)[names(tmpTable)=="V7"]="inf1"
    names(tmpTable)[names(tmpTable)=="V8"]="sup1"
    names(tmpTable)[names(tmpTable)=="V9"]="inf2"
    names(tmpTable)[names(tmpTable)=="V10"]="sup2"
  }
  # else : on suppose que les noms des colonnes sont bien

  if(!is.element("orient",names(tmpTable))){
    tmpTable$orient=ifelse(tmpTable$inf2<tmpTable$sup2,1,-1)
  }
  
  ## met les inf2<sup2 :
  inf2=pmin(tmpTable$inf2,tmpTable$sup2)
  sup2=pmax(tmpTable$inf2,tmpTable$sup2)
  tmpTable$inf2=inf2
  tmpTable$sup2=sup2

  removeList=unlist(by(tmpTable,couples,function(x) {
    ## pour chaque couple query-target
    l1=getRemoveList(x[x$orient==1,c("inf1","sup1","inf2","sup2","id")])
    l2=getRemoveList(x[x$orient==-1,c("inf1","sup1","inf2","sup2","id")])
    removeId=c(l1,l2)
    return(removeId)
  }))

  if(rename){
    return(tmpTable[!is.element(tmpTable$id,removeList),])
  }
  else{
    return(blastTable[!is.element(blastTable$id,removeList),])
  }
  
}


# fonction pour un tableau correspondant a un couple query-target et une seule orientation
getRemoveList=function(tab){

  removeId=NULL
  if(nrow(tab)>1){
    ## enlever les lignes identiques
    v=duplicated(tab[,c("inf1","sup1","inf2","sup2")])
    removeId=tab$id[v]
    ## pour les autres :
    autres=tab[!v,]
    if(nrow(autres)>1){
      for (i in 1:nrow(autres)){
        if(isEmbedded(autres[i,],autres[-i,])){
          removeId=c(removeId,autres$id[i])
        }
      }
    }
  }

  return(removeId)
}

isEmbedded=function(positions,restTable){
  e=sum(restTable$inf1<=positions$inf1 & restTable$sup1>=positions$sup1 & restTable$inf2<=positions$inf2 & restTable$sup2>=positions$sup2)
  if(e>0){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

getOverlapList=function(tab,fun="overlapsOn1",seuil=0){

  overlapId=NULL
  if(nrow(tab)>1){
    for (i in 1:nrow(tab)){
      if(eval(call(fun,tab[i,],tab[-i,],seuil))){
        overlapId=c(overlapId,tab$id[i])
      }
    }
  }

  return(overlapId)
}

overlapsOn1=function(positions,restTable,mini=0){
  overSize=pmin(restTable$sup1,positions$sup1)-pmax(restTable$inf1,positions$inf1)+1
    if(sum(overSize>mini)>0){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
}

overlapsOn1ButNot2=function(positions,restTable,mini=0){
  overSize1=pmin(restTable$sup1,positions$sup1)-pmax(restTable$inf1,positions$inf1)+1
  overSize2=pmin(restTable$sup2,positions$sup2)-pmax(restTable$inf2,positions$inf2)+1
  if(sum(overSize1>mini & overSize2<0)>0){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}



getWellOrderedPairs=function(tab,orient=1,seuil=0){
# pairs of hits with overlap1>seuil (but not one embedded in the other), no overlap on 2, and reverse order between 1 and 2 if orient=1 :
  # if orient=1 : the first hit of the pair is located before the second on the read, but after on the scaffolfd
  # if orient=-1 : same order

  # returns a table of pairs of hits
  # warning : can be be big when combining all pairs of overlapping reads

  res=by(tab,tab$couple,function(x) {
    if(nrow(x)>1){
      id1=NULL;id2=NULL
      t1=x[order(x$inf1,-x$sup1),c("inf1","sup1","inf2","sup2","id")]
      for (i in 1:(nrow(t1)-1)){
        sup1=t1$sup1[i]
        t2=t1[(i+1):nrow(t1),]
        t2=t2[t2$inf1<=sup1-seuil & t2$sup1>sup1,] ## overlap sur 1 d'au moins seuil et pas embedded
        if(nrow(t2)>0){
          for (j in 1:nrow(t2)){
            overSize2=min(t2$sup2[j],t1$sup2[i])-max(t2$inf2[j],t1$inf2[i])+1
            if((overSize2<0) & ((orient==1 & t2$inf2[j]<t1$inf2[i]) | (orient==-1 & t2$inf2[j]>t1$inf2[i]))){
              ## bonne paire
              id1=c(id1,t1$id[i])
              id2=c(id2,t2$id[j])
            }
          }
        }
      }
      return(data.frame(id1,id2))
    }
  })

  id1=unlist(lapply(res,function(x) {x$id1}),use.names=F)
  id2=unlist(lapply(res,function(x) {x$id2}),use.names=F)
  idPairs=data.frame(id1,id2)

  finalTab=merge(idPairs,tab,by.x="id1",by.y="id")
  finalTab=merge(finalTab,tab[,c("id","pcid","length","inf1","sup1","inf2","sup2","evalue")],by.x="id2",by.y="id",suffixes=c(".1",".2"))

  finalTab$totalSize=pmax(finalTab$sup2.1,finalTab$sup2.2)-pmin(finalTab$inf2.1,finalTab$inf2.2)+1
  finalTab$drjSize=finalTab$sup1.1-finalTab$inf1.2
  
  return(finalTab)

}
